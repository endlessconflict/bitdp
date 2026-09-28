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
            vp: []u64,

            pub fn init(gpa: std.mem.Allocator, pattern: []const u8) !Aligner {
                std.debug.assert(pattern.len >= 1);
                // One spare bit: the bottom row's state is read at bit m.
                const words = pattern.len / 64 + 1;
                const cpl = try gpa.alloc(u64, ncp * 256 * words);
                errdefer gpa.free(cpl);
                fillPlanes(pattern, words, cpl);
                return .{ .m = pattern.len, .words = words, .cpl = cpl, .vp = try gpa.alloc(u64, nl * words) };
            }

            pub fn deinit(a: *Aligner, gpa: std.mem.Allocator) void {
                gpa.free(a.cpl);
                gpa.free(a.vp);
            }

            pub fn distance(a: *Aligner, text: []const u8) i64 {
                return columns(null, a.words, a.m, a.cpl, a.vp, text);
            }
        };

        /// Pattern length 1..63, one machine word, no allocation.
        pub fn distance(pattern: []const u8, text: []const u8) i64 {
            std.debug.assert(pattern.len >= 1 and pattern.len <= 63);
            var cpl: [ncp * 256]u64 = undefined;
            var vp: [nl]u64 = undefined;
            fillPlanes(pattern, 1, &cpl);
            return columns(1, 1, pattern.len, &cpl, &vp, text);
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

        /// The DP over text columns. Within a column the derived program runs
        /// word by word, low rows first; additions pass their carry and
        /// shifts their top bit on to the next word.
        fn columns(comptime ct_words: ?usize, rt_words: usize, m: usize, cpl: []const u64, vp: []u64, text: []const u8) i64 {
            @setEvalBranchQuota(1 << 20);
            const words = ct_words orelse rt_words;
            for (0..nl) |t| {
                const bit = scheme.gap >= plan.vals[t + 1];
                @memset(vp[t * words ..][0..words], if (bit != plan.inp_pol[t]) ~@as(u64, 0) else 0);
            }
            const score_word = m / 64;
            const shift: u6 = @intCast(m % 64);
            var score: i64 = @as(i64, @intCast(m)) * scheme.gap;
            for (text) |c| {
                var carry: [plan.len]u64 = undefined;
                var top: [plan.len]u64 = undefined;
                for (0..words) |w| {
                    var r: [plan.len]u64 = undefined;
                    inline for (plan.nodes[0..plan.len], 0..) |n, i| {
                        r[i] = switch (n.op) {
                            .input => if (n.a < nl) vp[n.a * words + w] else cpl[((n.a - nl) * 256 + c) * words + w],
                            .zero => 0,
                            .ones => ~@as(u64, 0),
                            .not => ~r[n.a],
                            .@"and" => r[n.a] & r[n.b],
                            .@"or" => r[n.a] | r[n.b],
                            .xor => r[n.a] ^ r[n.b],
                            .add, .add1 => blk: {
                                const cin: u64 = if (w == 0) @intFromBool(n.op == .add1) else carry[i];
                                const s1 = @addWithOverflow(r[n.a], r[n.b]);
                                const s2 = @addWithOverflow(s1[0], cin);
                                carry[i] = s1[1] | s2[1];
                                break :blk s2[0];
                            },
                            .shl0, .shl1 => blk: {
                                const fill: u64 = if (w == 0) @intFromBool(n.op == .shl1) else top[i];
                                top[i] = r[n.a] >> 63;
                                break :blk (r[n.a] << 1) | fill;
                            },
                        };
                    }
                    if (w == score_word) {
                        var idx: usize = 0;
                        inline for (0..nl) |t| {
                            const bit = (r[plan.out_u[t]] >> shift) & 1;
                            idx += bit ^ @intFromBool(plan.lvl_pol[t + 1]);
                        }
                        score += plan.vals[idx];
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
