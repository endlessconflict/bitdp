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
pub const deriveWith = derive_mod.deriveWith;
pub const Options = derive_mod.Options;
pub const derive_opcost = derive_mod.opCost;

/// Alignment cost kernel for `scheme` (global or search mode).
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
                // One spare row: the bottom row's state is read at row m.
                const words = pattern.len / rows + 1;
                const cpl = try gpa.alloc(u64, plane_words * words);
                errdefer gpa.free(cpl);
                fillPlanes(pattern, words, cpl);
                return .{ .m = pattern.len, .words = words, .cpl = cpl, .vp = try gpa.alloc(@Vector(1, u64), nl * words) };
            }

            pub fn deinit(a: *Aligner, gpa: std.mem.Allocator) void {
                gpa.free(a.cpl);
                gpa.free(a.vp);
            }

            pub fn distance(a: *Aligner, text: []const u8) i64 {
                return columns(1, false, null, a.words, .{a.m}, .{a.cpl}, &.{}, a.vp, .{text}, null)[0];
            }
        };

        /// Pattern length 1..62, one machine word, no allocation.
        pub fn distance(pattern: []const u8, text: []const u8) i64 {
            std.debug.assert(pattern.len >= 1 and pattern.len < rows);
            var cpl: [plane_words]u64 = undefined;
            var vp: [nl]@Vector(1, u64) = undefined;
            fillPlanes(pattern, 1, &cpl);
            return columns(1, false, 1, 1, .{pattern.len}, .{&cpl}, &.{}, &vp, .{text}, null)[0];
        }

        /// SIMD lanes: one independent alignment per lane.
        pub const lanes = std.simd.suggestVectorLength(u64) orelse 4;
        const V = @Vector(lanes, u64);

        /// Up to `lanes` patterns prepared once, aligned together against a
        /// shared text: every lane reads the same text character, so each
        /// input plane is one vector load.
        pub const Group = struct {
            words: usize,
            ms: [lanes]usize,
            planes: []V,
            vp: []V,

            pub fn init(gpa: std.mem.Allocator, patterns: []const []const u8) !Group {
                std.debug.assert(patterns.len >= 1 and patterns.len <= lanes);
                var g: Group = .{ .words = 1, .ms = undefined, .planes = undefined, .vp = undefined };
                for (0..lanes) |l| {
                    g.ms[l] = patterns[@min(l, patterns.len - 1)].len; // pad with the last pattern
                    g.words = @max(g.words, g.ms[l] / rows + 1);
                }
                g.planes = try gpa.alloc(V, plane_words * g.words);
                errdefer gpa.free(g.planes);
                const one = try gpa.alloc(u64, plane_words * g.words);
                defer gpa.free(one);
                inline for (0..lanes) |l| {
                    fillPlanes(patterns[@min(l, patterns.len - 1)], g.words, one);
                    for (g.planes, one) |*v, x| v[l] = x;
                }
                g.vp = try gpa.alloc(V, nl * g.words);
                return g;
            }

            pub fn deinit(g: *Group, gpa: std.mem.Allocator) void {
                gpa.free(g.planes);
                gpa.free(g.vp);
            }

            /// Search mode: every end position in `text` where some pattern of
            /// the group aligns with cost at most `max_cost`, appended to `hits`.
            pub fn scan(g: *Group, gpa: std.mem.Allocator, text: []const u8, max_cost: i64, hits: *std.ArrayList(Hit)) !void {
                comptime std.debug.assert(scheme.mode == .search);
                var sink: Sink = .{ .gpa = gpa, .hits = hits, .max_cost = max_cost };
                _ = columns(lanes, true, null, g.words, g.ms, undefined, g.planes, g.vp, @splat(text), &sink);
                if (sink.err) |e| return e;
            }

            /// Cost of every pattern in the group against `text`.
            pub fn distances(g: *Group, text: []const u8) [lanes]i64 {
                return columns(lanes, true, null, g.words, g.ms, undefined, g.planes, g.vp, @splat(text), null);
            }
        };

        /// A search-mode match: pattern `lane` of a group ends at text position
        /// `end` (exclusive) with alignment cost `cost`.
        pub const Hit = struct { lane: u8, end: usize, cost: i64 };

        const Sink = struct {
            gpa: std.mem.Allocator,
            hits: *std.ArrayList(Hit),
            max_cost: i64,
            err: ?anyerror = null,
        };

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
                    words = @max(words, ms[l] / rows + 1);
                }
                const buf = try gpa.alloc(u64, lanes * plane_words * words);
                defer gpa.free(buf);
                const vp = try gpa.alloc(V, nl * words);
                defer gpa.free(vp);
                var cpls: [lanes][]const u64 = undefined;
                for (0..lanes) |l| {
                    const mine = buf[l * plane_words * words ..][0 .. plane_words * words];
                    fillPlanes(patterns[@min(start + l, patterns.len - 1)], words, mine);
                    cpls[l] = mine;
                }
                const r = columns(lanes, false, null, words, ms, cpls, &.{}, vp, ts, null);
                for (0..@min(lanes, patterns.len - start)) |l| out[start + l] = r[l];
            }
        }

        /// Pattern rows per 64-bit word. Bit 63 is left as a carry catcher: the
        /// carry out of rows 0..62 is ((a + b) ^ a ^ b) >> 63 whatever bit 63
        /// holds, which costs three instructions where a full 64-bit carry-out
        /// would need an emulated unsigned compare on most SIMD units.
        const rows = 63;

        /// Text bytes map to slots: one per alphabet character, plus one for
        /// every other byte, which counts as the costliest substitution.
        const nslot = scheme.alphabet.len + 1;
        const slot_of: [256]u8 = blk: {
            var t = [_]u8{scheme.alphabet.len} ** 256;
            for (scheme.alphabet, 0..) |c, i| t[c] = i;
            break :blk t;
        };
        const plane_words = ncp * nslot;

        /// Class planes of one pattern: word w of plane j for text slot s is at
        /// (j * nslot + s) * words + w, with bit i set when the cost of pattern[i]
        /// against that slot's character is at least cls.c[j+1].
        fn fillPlanes(pattern: []const u8, words: usize, cpl: []u64) void {
            @memset(cpl, 0);
            for (0..nslot) |s| for (pattern, 0..) |c, i| {
                const cost = if (s < scheme.alphabet.len) scheme.cost(c, scheme.alphabet[s]) else plan.cls.c[plan.cls.n - 1];
                inline for (0..ncp) |j| {
                    if (cost >= plan.cls.c[j + 1]) cpl[(j * nslot + s) * words + i / rows] |= @as(u64, 1) << @intCast(i % rows);
                }
            };
            inline for (0..ncp) |j| {
                if (plan.inp_pol[nl + j]) for (cpl[j * nslot * words ..][0 .. nslot * words]) |*x| {
                    x.* = ~x.*;
                };
            }
        }

        /// The DP over text columns, for L independent alignments in the L
        /// lanes of a vector. Within a column the derived program runs word by
        /// word, low rows first; additions pass their carry and shifts their
        /// top bit on to the next word. A lane stops scoring once its text ends.
        /// Class planes come per lane (`cpls`) or, when every lane reads the
        /// same text (`shared`), pre-interleaved (`gplanes`).
        fn columns(
            comptime L: usize,
            comptime shared: bool,
            comptime ct_words: ?usize,
            rt_words: usize,
            ms: [L]usize,
            cpls: [L][]const u64,
            gplanes: []const @Vector(L, u64),
            vp: []@Vector(L, u64),
            texts: [L][]const u8,
            sink: ?*Sink,
        ) [L]i64 {
            @setEvalBranchQuota(1 << 20);
            const W = @Vector(L, u64);
            const S = @Vector(L, i64);
            const words = ct_words orelse rt_words;
            const ones: W = @splat(~@as(u64, 0));
            const zeros: W = @splat(0);
            for (0..nl) |t| {
                const bit = scheme.gap >= plan.vals[t + 1];
                @memset(vp[t * words ..][0..words], if (bit != plan.inp_pol[t]) ones else zeros);
            }
            var score_word: W = undefined;
            var shift: @Vector(L, u6) = undefined;
            var score: S = undefined;
            var n_max: usize = 0;
            var sw_min: usize = std.math.maxInt(usize);
            var sw_max: usize = 0;
            inline for (0..L) |l| {
                score_word[l] = ms[l] / rows;
                shift[l] = @intCast(ms[l] % rows);
                sw_min = @min(sw_min, ms[l] / rows);
                sw_max = @max(sw_max, ms[l] / rows);
                score[l] = @as(i64, @intCast(ms[l])) * scheme.gap;
                n_max = @max(n_max, texts[l].len);
            }
            var best = score; // search mode: minimum over end positions, D(m, 0) included
            for (0..n_max) |j| {
                var cs: [L]usize = undefined;
                var active: @Vector(L, bool) = @splat(true);
                if (shared) {
                    cs[0] = slot_of[texts[0][j]];
                } else inline for (0..L) |l| {
                    active[l] = j < texts[l].len;
                    cs[l] = if (active[l]) slot_of[texts[l][j]] else 0;
                }
                var carry: [plan.len]W = undefined;
                var top: [plan.len]W = undefined;
                for (0..words) |w| {
                    var r: [plan.len]W = undefined;
                    inline for (plan.nodes[0..plan.len], 0..) |n, i| {
                        r[i] = switch (n.op) {
                            .input => if (n.a < nl) vp[n.a * words + w] else if (shared)
                                gplanes[((n.a - nl) * nslot + cs[0]) * words + w]
                            else blk: {
                                var x: W = undefined;
                                inline for (0..L) |l| x[l] = cpls[l][((n.a - nl) * nslot + cs[l]) * words + w];
                                break :blk x;
                            },
                            .zero => zeros,
                            .ones => ones,
                            .not => ~r[n.a],
                            .@"and" => r[n.a] & r[n.b],
                            .@"or" => r[n.a] | r[n.b],
                            .xor => r[n.a] ^ r[n.b],
                            .add, .add1 => blk: {
                                const cin: W = if (w == 0) @splat(@intFromBool(n.op == .add1)) else carry[i];
                                const sum = r[n.a] +% r[n.b] +% cin;
                                carry[i] = (sum ^ r[n.a] ^ r[n.b]) >> @splat(rows);
                                break :blk sum;
                            },
                            .shl0, .shl1 => blk: {
                                const fill: W = if (w == 0) @splat(@intFromBool(n.op == .shl1)) else top[i];
                                top[i] = (r[n.a] >> @splat(rows - 1)) & @as(W, @splat(1));
                                break :blk (r[n.a] << @splat(1)) | fill;
                            },
                        };
                    }
                    // Bottom-row state, read at row m of the lanes whose m falls in this word.
                    if (w >= sw_min and w <= sw_max) {
                        const here = score_word == @as(W, @splat(w));
                        var inc: S = @splat(plan.vals[0]);
                        inline for (0..nl) |t| {
                            const pol: W = @splat(@intFromBool(plan.lvl_pol[t + 1]));
                            const bit: S = @intCast(((r[plan.out_u[t]] >> shift) & @as(W, @splat(1))) ^ pol);
                            inc += bit * @as(S, @splat(plan.vals[t + 1] - plan.vals[t]));
                        }
                        score += @select(i64, @select(bool, here, active, @as(@Vector(L, bool), @splat(false))), inc, @as(S, @splat(0)));
                    }
                    inline for (0..nl) |t| vp[t * words + w] = r[plan.out_v[t]];
                }
                if (scheme.mode == .search) {
                    best = @select(i64, active, @min(best, score), best);
                    if (sink) |k| {
                        const hit = @select(bool, active, score <= @as(S, @splat(k.max_cost)), @as(@Vector(L, bool), @splat(false)));
                        if (@reduce(.Or, hit)) inline for (0..L) |l| {
                            if (hit[l]) k.hits.append(k.gpa, .{ .lane = l, .end = j + 1, .cost = score[l] }) catch |e| {
                                k.err = e;
                            };
                        };
                    }
                }
            }
            return if (scheme.mode == .search) best else score;
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
        const m = rnd.intRangeAtMost(usize, 1, 62);
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

test "groups: several patterns against one shared text match the scalar oracle" {
    const gpa = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(8);
    const rnd = prng.random();
    var buf: [301]i64 = undefined;
    inline for ([_]Scheme{ .{}, schemes.ts_tv }) |s| {
        const K = Kernel(s);
        var store: [K.lanes][300]u8 = undefined;
        var ps: [K.lanes][]const u8 = undefined;
        for (0..20) |_| {
            const count = rnd.intRangeAtMost(usize, 1, K.lanes);
            for (0..count) |i| {
                const m = rnd.intRangeAtMost(usize, 1, 300);
                for (store[i][0..m]) |*x| x.* = "ACGTN"[rnd.intRangeLessThan(usize, 0, 5)];
                ps[i] = store[i][0..m];
            }
            var text: [300]u8 = undefined;
            const n = rnd.intRangeAtMost(usize, 0, 300);
            for (text[0..n]) |*x| x.* = "ACGT"[rnd.int(u2)];
            var g = try K.Group.init(gpa, ps[0..count]);
            defer g.deinit(gpa);
            const r = g.distances(text[0..n]);
            for (0..count) |i| try std.testing.expectEqual(reference.scalar(s, ps[i], text[0..n], &buf), r[i]);
        }
    }
}

test "search mode (pattern against the best substring) matches the scalar oracle" {
    const gpa = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(12);
    const rnd = prng.random();
    var buf: [201]i64 = undefined;
    inline for ([_]Scheme{
        .{ .mode = .search },
        .{ .match = -2, .mismatch = 3, .gap = 5, .mode = .search },
        .{ .sub = &schemes.tsTvCost, .gap = 2, .mode = .search },
    }) |s| {
        const K = Kernel(s);
        var p: [200]u8 = undefined;
        var t: [300]u8 = undefined;
        for (0..400) |_| {
            const m = rnd.intRangeAtMost(usize, 1, 200);
            const n = rnd.intRangeAtMost(usize, 0, 300);
            for (p[0..m]) |*x| x.* = "ACGT"[rnd.int(u2)];
            for (t[0..n]) |*x| x.* = "ACGT"[rnd.int(u2)];
            // Plant the pattern (mutated) inside the text half of the time.
            if (n > m and rnd.boolean()) {
                const at = rnd.uintLessThan(usize, n - m);
                @memcpy(t[at..][0..m], p[0..m]);
                t[at + m / 2] = 'A';
            }
            const want = reference.scalar(s, p[0..m], t[0..n], &buf);
            var a = try K.Aligner.init(gpa, p[0..m]);
            defer a.deinit(gpa);
            try std.testing.expectEqual(want, a.distance(t[0..n]));
            if (m < 63) try std.testing.expectEqual(want, K.distance(p[0..m], t[0..n]));
        }
    }
}

test "group scan reports exactly the end positions the scalar DP finds" {
    const gpa = std.testing.allocator;
    const s: Scheme = .{ .sub = &schemes.tsTvCost, .gap = 2, .mode = .search };
    const K = Kernel(s);
    var prng = std.Random.DefaultPrng.init(13);
    const rnd = prng.random();
    var text: [2000]u8 = undefined;
    for (&text) |*x| x.* = "ACGT"[rnd.int(u2)];
    var store: [K.lanes][20]u8 = undefined;
    var ps: [K.lanes][]const u8 = undefined;
    for (0..K.lanes) |l| {
        const at = rnd.uintLessThan(usize, text.len - 20);
        @memcpy(&store[l], text[at..][0..20]);
        store[l][rnd.uintLessThan(usize, 20)] = 'T';
        ps[l] = &store[l];
    }
    var g = try K.Group.init(gpa, &ps);
    defer g.deinit(gpa);
    var hits: std.ArrayList(K.Hit) = .empty;
    defer hits.deinit(gpa);
    try g.scan(gpa, &text, 4, &hits);
    // Scalar reference: the last DP row, column by column, is the best cost of
    // the pattern ending at each text position.
    var expected: usize = 0;
    for (0..K.lanes) |l| {
        var col: [21]i64 = undefined;
        for (&col, 0..) |*x, i| x.* = @as(i64, @intCast(i)) * s.gap;
        for (text, 1..) |c, e| {
            var diag = col[0];
            col[0] = 0;
            for (ps[l], 1..) |pc, i| {
                const cell = @min(diag + s.cost(pc, c), @min(col[i] + s.gap, col[i - 1] + s.gap));
                diag = col[i];
                col[i] = cell;
            }
            if (col[20] <= 4) {
                expected += 1;
                var found = false;
                for (hits.items) |h| found = found or (h.lane == l and h.end == e and h.cost == col[20]);
                try std.testing.expect(found);
            }
        }
    }
    try std.testing.expectEqual(expected, hits.items.len);
}
