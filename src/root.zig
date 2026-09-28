//! bitdp: bit-parallel alignment kernels derived at compile time from the
//! scoring scheme. See derive.zig for how the derivation works.

const std = @import("std");
const derive_mod = @import("derive.zig");
pub const reference = @import("reference.zig");

pub const Scheme = derive_mod.Scheme;
pub const Plan = derive_mod.Plan;
pub const Op = derive_mod.Op;
pub const derive = derive_mod.derive;

/// Global alignment cost kernel for `scheme`, pattern length 1..64.
pub fn Kernel(comptime scheme: Scheme) type {
    const plan = comptime derive(scheme);
    const nl = plan.k - 1;
    return struct {
        pub const ops = plan.cost;
        pub const chains = plan.chains;
        pub const levels = nl;

        pub fn distance(pattern: []const u8, text: []const u8) i64 {
            const m = pattern.len;
            std.debug.assert(m >= 1 and m <= 64);
            var peq = [_]u64{0} ** 256;
            for (pattern, 0..) |c, i| peq[c] |= @as(u64, 1) << @intCast(i);

            // Column 0: every Delta-v equals gap.
            var vp: [nl]u64 = undefined;
            inline for (0..nl) |t| {
                const bit = scheme.gap >= plan.vals[t + 1];
                vp[t] = if (bit != plan.in_pol) ~@as(u64, 0) else 0;
            }
            const shift: u6 = @intCast(m - 1);
            var score: i64 = @as(i64, @intCast(m)) * scheme.gap;

            for (text) |c| {
                var r: [plan.len]u64 = undefined;
                inline for (plan.nodes[0..plan.len], 0..) |n, i| {
                    r[i] = switch (n.op) {
                        .input => if (n.a < nl) vp[n.a] else peq[c],
                        .zero => 0,
                        .ones => ~@as(u64, 0),
                        .not => ~r[n.a],
                        .@"and" => r[n.a] & r[n.b],
                        .@"or" => r[n.a] | r[n.b],
                        .xor => r[n.a] ^ r[n.b],
                        .add => r[n.a] +% r[n.b],
                        .add1 => r[n.a] +% r[n.b] +% 1,
                        .shl0 => r[n.a] << 1,
                        .shl1 => (r[n.a] << 1) | 1,
                    };
                }
                var idx: usize = 0;
                inline for (0..nl) |t| {
                    const bit = (r[plan.out_u[t]] >> shift) & 1;
                    idx += bit ^ @intFromBool(plan.lvl_pol[t + 1]);
                }
                score += plan.vals[idx];
                inline for (0..nl) |t| vp[t] = r[plan.out_v[t]];
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
        const m = rnd.intRangeAtMost(usize, 1, 64);
        const n = rnd.intRangeAtMost(usize, 0, t.len);
        for (p[0..m]) |*x| x.* = "ACGT"[rnd.int(u2)];
        for (t[0..n]) |*x| x.* = "ACGT"[rnd.int(u2)];
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

test "MVE gate: derived unit-cost kernel within 1.5x of Myers' op count" {
    const K = Kernel(.{});
    try std.testing.expect(K.chains == 1);
    try std.testing.expect(2 * K.ops <= 3 * reference.myers_ops);
}
