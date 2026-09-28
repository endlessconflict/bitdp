//! bitdp: bit-parallel alignment kernels derived at compile time from the
//! scoring scheme. See derive.zig for how the derivation works.

const std = @import("std");
const derive_mod = @import("derive.zig");
pub const reference = @import("reference.zig");
pub const schemes = @import("schemes.zig");

pub const Scheme = derive_mod.Scheme;
pub const Plan = derive_mod.Plan;
pub const Op = derive_mod.Op;
pub const derive = derive_mod.derive;
pub const derive_opcost = derive_mod.opCost;

/// Global alignment cost kernel for `scheme`.
pub fn Kernel(comptime scheme: Scheme) type {
    const plan = comptime derive(scheme);
    const nl = plan.k - 1;
    const ncp = plan.cls.n - 1;
    return struct {
        pub const ops = plan.cost;
        pub const chains = plan.chains;
        pub const levels = nl;

        /// Pattern prepared once for texts of any length; any pattern length.
        pub const Aligner = struct {
            m: usize,
            words: usize,
            cpl: []u64,
            vp: []@Vector(1, u64),

            pub fn init(gpa: std.mem.Allocator, pattern: []const u8) !Aligner {
                std.debug.assert(pattern.len >= 1);
                // One spare bit: the bottom row's state is read at bit m.
                const words = pattern.len / 64 + 1;
                const cpl = try gpa.alloc(u64, ncp * 256 * words);
                errdefer gpa.free(cpl);
                fillPlanes(pattern, words, cpl);
                return .{ .m = pattern.len, .words = words, .cpl = cpl, .vp = try gpa.alloc(@Vector(1, u64), nl * words) };
            }

            pub fn deinit(a: *Aligner, gpa: std.mem.Allocator) void {
                gpa.free(a.cpl);
                gpa.free(a.vp);
            }

            pub fn distance(a: *Aligner, text: []const u8) i64 {
                return columns(1, null, a.words, .{a.m}, .{a.cpl}, a.vp, .{text})[0];
            }
        };

        /// Pattern length 1..63, one machine word, no allocation.
        pub fn distance(pattern: []const u8, text: []const u8) i64 {
            std.debug.assert(pattern.len >= 1 and pattern.len <= 63);
            var cpl: [ncp * 256]u64 = undefined;
            var vp: [nl]@Vector(1, u64) = undefined;
            fillPlanes(pattern, 1, &cpl);
            return columns(1, 1, 1, .{pattern.len}, .{&cpl}, &vp, .{text})[0];
        }

        /// SIMD lanes used by `distances`: one independent alignment per lane.
        pub const lanes = std.simd.suggestVectorLength(u64) orelse 4;

        /// Costs of many (pattern, text) pairs, `lanes` pairs at a time: each
        /// SIMD lane runs the same derived program on its own alignment.
        /// Pairs with similar pattern lengths batch best.
        pub fn distances(gpa: std.mem.Allocator, patterns: []const []const u8, texts: []const []const u8, out: []i64) !void {
            std.debug.assert(patterns.len == texts.len and out.len == patterns.len);
            var start: usize = 0;
            while (start < patterns.len) : (start += lanes) {
                var ms: [lanes]usize = undefined;
                var ts: [lanes][]const u8 = undefined;
                var words: usize = 1;
                for (0..lanes) |l| {
                    const i = @min(start + l, patterns.len - 1); // pad the tail with the last pair
                    ms[l] = patterns[i].len;
                    ts[l] = texts[i];
                    words = @max(words, ms[l] / 64 + 1);
                }
                const buf = try gpa.alloc(u64, lanes * ncp * 256 * words);
                defer gpa.free(buf);
                const vp = try gpa.alloc(@Vector(lanes, u64), nl * words);
                defer gpa.free(vp);
                var cpls: [lanes][]const u64 = undefined;
                for (0..lanes) |l| {
                    const slot = buf[l * ncp * 256 * words ..][0 .. ncp * 256 * words];
                    fillPlanes(patterns[@min(start + l, patterns.len - 1)], words, slot);
                    cpls[l] = slot;
                }
                const r = columns(lanes, null, words, ms, cpls, vp, ts);
                for (0..@min(lanes, patterns.len - start)) |l| out[start + l] = r[l];
            }
        }

        /// Class planes: word w of plane j for text byte c is at
        /// (j * 256 + c) * words + w; bit i set when cost(pattern[i], c) >= cls.c[j+1].
        fn fillPlanes(pattern: []const u8, words: usize, cpl: []u64) void {
            @memset(cpl, 0);
            if (ncp == 0) return;
            const bit = struct {
                fn f(i: usize) u64 {
                    return @as(u64, 1) << @intCast(i % 64);
                }
            }.f;
            if (scheme.sub == null) {
                for (pattern, 0..) |c, i| cpl[@as(usize, c) * words + i / 64] |= bit(i);
                // That marked equal bytes. [cost >= c1] means "mismatch" when mismatch > match.
                if (scheme.mismatch > scheme.match) for (cpl) |*x| {
                    x.* = ~x.*;
                };
            } else {
                for (scheme.alphabet) |a| for (pattern, 0..) |c, i| {
                    const cost = scheme.cost(c, a);
                    inline for (0..ncp) |j| {
                        if (cost >= plan.cls.c[j + 1]) cpl[(j * 256 + a) * words + i / 64] |= bit(i);
                    }
                };
            }
            inline for (0..ncp) |j| {
                if (plan.inp_pol[nl + j]) for (cpl[j * 256 * words ..][0 .. 256 * words]) |*x| {
                    x.* = ~x.*;
                };
            }
        }

        /// The DP over text columns, for L independent alignments in the L
        /// lanes of a vector. Within a column the derived program runs word by
        /// word, low rows first; additions pass their carry and shifts their
        /// top bit on to the next word. A lane stops scoring once its text ends.
        fn columns(
            comptime L: usize,
            comptime ct_words: ?usize,
            rt_words: usize,
            ms: [L]usize,
            cpls: [L][]const u64,
            vp: []@Vector(L, u64),
            texts: [L][]const u8,
        ) [L]i64 {
            @setEvalBranchQuota(1 << 20);
            const V = @Vector(L, u64);
            const S = @Vector(L, i64);
            const words = ct_words orelse rt_words;
            const ones: V = @splat(~@as(u64, 0));
            const zeros: V = @splat(0);
            for (0..nl) |t| {
                const bit = scheme.gap >= plan.vals[t + 1];
                @memset(vp[t * words ..][0..words], if (bit != plan.inp_pol[t]) ones else zeros);
            }
            var score_word: V = undefined;
            var shift: @Vector(L, u6) = undefined;
            var score: S = undefined;
            var n_max: usize = 0;
            inline for (0..L) |l| {
                score_word[l] = ms[l] / 64;
                shift[l] = @intCast(ms[l] % 64);
                score[l] = @as(i64, @intCast(ms[l])) * scheme.gap;
                n_max = @max(n_max, texts[l].len);
            }
            for (0..n_max) |j| {
                var cs: [L]usize = undefined;
                var active: @Vector(L, bool) = undefined;
                inline for (0..L) |l| {
                    active[l] = j < texts[l].len;
                    cs[l] = if (active[l]) texts[l][j] else 0;
                }
                var carry: [plan.len]V = undefined;
                var top: [plan.len]V = undefined;
                for (0..words) |w| {
                    var r: [plan.len]V = undefined;
                    inline for (plan.nodes[0..plan.len], 0..) |n, i| {
                        r[i] = switch (n.op) {
                            .input => if (n.a < nl) vp[n.a * words + w] else blk: {
                                var x: V = undefined;
                                inline for (0..L) |l| x[l] = cpls[l][((n.a - nl) * 256 + cs[l]) * words + w];
                                break :blk x;
                            },
                            .zero => zeros,
                            .ones => ones,
                            .not => ~r[n.a],
                            .@"and" => r[n.a] & r[n.b],
                            .@"or" => r[n.a] | r[n.b],
                            .xor => r[n.a] ^ r[n.b],
                            .add, .add1 => blk: {
                                const cin: V = if (w == 0) @splat(@intFromBool(n.op == .add1)) else carry[i];
                                const s1 = @addWithOverflow(r[n.a], r[n.b]);
                                const s2 = @addWithOverflow(s1[0], cin);
                                const c1: V = s1[1];
                                const c2: V = s2[1];
                                carry[i] = c1 | c2;
                                break :blk s2[0];
                            },
                            .shl0, .shl1 => blk: {
                                const fill: V = if (w == 0) @splat(@intFromBool(n.op == .shl1)) else top[i];
                                top[i] = r[n.a] >> @splat(63);
                                break :blk (r[n.a] << @splat(1)) | fill;
                            },
                        };
                    }
                    // Bottom-row state, read at bit m of the lanes whose m falls in this word.
                    const here = score_word == @as(V, @splat(w));
                    if (@reduce(.Or, here)) {
                        var inc: S = @splat(plan.vals[0]);
                        inline for (0..nl) |t| {
                            const pol: V = @splat(@intFromBool(plan.lvl_pol[t + 1]));
                            const bit: S = @intCast(((r[plan.out_u[t]] >> shift) & @as(V, @splat(1))) ^ pol);
                            inc += bit * @as(S, @splat(plan.vals[t + 1] - plan.vals[t]));
                        }
                        score += @select(i64, @select(bool, here, active, @as(@Vector(L, bool), @splat(false))), inc, @as(S, @splat(0)));
                    }
                    inline for (0..nl) |t| vp[t * words + w] = r[plan.out_v[t]];
                }
            }
            return score;
        }
    };
}

// ---------------------------------------------------------------- tests

fn checkRandom(comptime s: Scheme, pairs: usize, seed: u64) !void {
    const K = Kernel(s);
    var prng = std.Random.DefaultPrng.init(seed);
    const rnd = prng.random();
    var buf: [65]i64 = undefined;
    var p: [64]u8 = undefined;
    var t: [200]u8 = undefined;
    for (0..pairs) |_| {
        const m = rnd.intRangeAtMost(usize, 1, 63);
        const n = rnd.intRangeAtMost(usize, 0, t.len);
        for (p[0..m]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
        for (t[0..n]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
        try std.testing.expectEqual(
            reference.scalar(s, p[0..m], t[0..n], &buf),
            K.distance(p[0..m], t[0..n]),
        );
    }
}

test "unit edit distance: random pairs match the scalar oracle" {
    try checkRandom(.{}, 20_000, 1);
}

test "unit edit distance: exhaustive over {A,C,G}, lengths 1..4 x 0..4" {
    const K = Kernel(.{});
    var buf: [8]i64 = undefined;
    var pat: [4]u8 = undefined;
    var txt: [4]u8 = undefined;
    for (1..5) |m| for (0..5) |n| {
        var pc: usize = 0;
        while (pc < std.math.pow(usize, 3, m)) : (pc += 1) {
            var x = pc;
            for (pat[0..m]) |*ch| {
                ch.* = "ACG"[x % 3];
                x /= 3;
            }
            var tc: usize = 0;
            while (tc < std.math.pow(usize, 3, n)) : (tc += 1) {
                var y = tc;
                for (txt[0..n]) |*ch| {
                    ch.* = "ACG"[y % 3];
                    y /= 3;
                }
                try std.testing.expectEqual(
                    reference.scalar(.{}, pat[0..m], txt[0..n], &buf),
                    K.distance(pat[0..m], txt[0..n]),
                );
            }
        }
    };
}

test "Myers baseline matches the scalar oracle" {
    var prng = std.Random.DefaultPrng.init(2);
    const rnd = prng.random();
    var buf: [65]i64 = undefined;
    var p: [64]u8 = undefined;
    var t: [100]u8 = undefined;
    for (0..5_000) |_| {
        const m = rnd.intRangeAtMost(usize, 1, 64);
        const n = rnd.intRangeAtMost(usize, 0, t.len);
        for (p[0..m]) |*x| x.* = "ACGT"[rnd.int(u2)];
        for (t[0..n]) |*x| x.* = "ACGT"[rnd.int(u2)];
        try std.testing.expectEqual(reference.scalar(.{}, p[0..m], t[0..n], &buf), reference.myers(p[0..m], t[0..n]));
    }
}

test "derived unit-cost kernel: one carry chain, no more word ops than Myers' hand-derived kernel" {
    const K = Kernel(.{});
    try std.testing.expect(K.chains == 1);
    try std.testing.expect(K.ops <= reference.myers_ops);
}

test "BitPAl weights (2, -3, -5) as costs: exact, and cheaper than BitPAl's 265 ops/word" {
    const s: Scheme = .{ .match = -2, .mismatch = 3, .gap = 5 };
    try checkRandom(s, 3_000, 3);
    // 265 is the operation count Loving et al. (2014) report for these weights.
    try std.testing.expect(Kernel(s).ops < 265);
}

test "three cost classes (transition 1, transversion 2): matches the scalar oracle" {
    try checkRandom(.{ .sub = &schemes.tsTvCost, .gap = 2 }, 5_000, 4);
}

test "multi-word patterns (up to 300) match the scalar oracle" {
    const gpa = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(5);
    const rnd = prng.random();
    var buf: [301]i64 = undefined;
    var p: [300]u8 = undefined;
    var t: [300]u8 = undefined;
    inline for ([_]Scheme{ .{}, .{ .match = -2, .mismatch = 3, .gap = 5 }, .{ .sub = &schemes.tsTvCost, .gap = 2 } }) |s| {
        for (0..300) |_| {
            const m = rnd.intRangeAtMost(usize, 1, 300);
            const n = rnd.intRangeAtMost(usize, 0, 300);
            for (p[0..m]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
            for (t[0..n]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
            var a = try Kernel(s).Aligner.init(gpa, p[0..m]);
            defer a.deinit(gpa);
            try std.testing.expectEqual(reference.scalar(s, p[0..m], t[0..n], &buf), a.distance(t[0..n]));
        }
    }
}

test "batched SIMD lanes match the scalar oracle (mixed lengths)" {
    const gpa = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(6);
    const rnd = prng.random();
    var buf: [201]i64 = undefined;
    inline for ([_]Scheme{ .{}, .{ .match = -2, .mismatch = 3, .gap = 5 } }) |s| {
        const K = Kernel(s);
        const count = 3 * K.lanes + 1;
        var store: [count][2][200]u8 = undefined;
        var ps: [count][]const u8 = undefined;
        var ts: [count][]const u8 = undefined;
        for (0..count) |i| {
            const m = rnd.intRangeAtMost(usize, 1, 200);
            const n = rnd.intRangeAtMost(usize, 0, 200);
            for (store[i][0][0..m]) |*x| x.* = "ACGT"[rnd.int(u2)];
            for (store[i][1][0..n]) |*x| x.* = "ACGT"[rnd.int(u2)];
            ps[i] = store[i][0][0..m];
            ts[i] = store[i][1][0..n];
        }
        var out: [count]i64 = undefined;
        try K.distances(gpa, &ps, &ts, &out);
        for (0..count) |i| try std.testing.expectEqual(reference.scalar(s, ps[i], ts[i], &buf), out[i]);
    }
}
