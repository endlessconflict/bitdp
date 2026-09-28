//! Reference implementations used as test oracles and benchmark baselines.

const std = @import("std");
const Scheme = @import("derive.zig").Scheme;

/// Plain O(mn) global alignment cost. The oracle every kernel is checked against.
pub fn scalar(s: Scheme, pattern: []const u8, text: []const u8, buf: []i64) i64 {
    const m = pattern.len;
    const col = buf[0 .. m + 1];
    for (col, 0..) |*x, i| x.* = @as(i64, @intCast(i)) * s.gap;
    for (text, 1..) |c, j| {
        var diag = col[0];
        col[0] = @as(i64, @intCast(j)) * s.gap;
        for (pattern, 1..) |pc, i| {
            const sub: i64 = s.cost(pc, c);
            const best = @min(diag + sub, @min(col[i] + s.gap, col[i - 1] + s.gap));
            diag = col[i];
            col[i] = best;
        }
    }
    return col[m];
}

/// Myers' bit-vector algorithm (myers1999bitvector), unit costs, global mode,
/// pattern length <= 64. Written from the recurrence as the hand-derived baseline.
/// 15 word operations per column (score update excluded).
pub fn myers(pattern: []const u8, text: []const u8) i64 {
    const m = pattern.len;
    std.debug.assert(m >= 1 and m <= 64);
    var peq = [_]u64{0} ** 256;
    for (pattern, 0..) |c, i| peq[c] |= @as(u64, 1) << @intCast(i);
    const top = @as(u64, 1) << @intCast(m - 1);
    var pv: u64 = ~@as(u64, 0);
    var mv: u64 = 0;
    var score: i64 = @intCast(m);
    for (text) |c| {
        const eq = peq[c];
        const xv = eq | mv;
        const xh = (((eq & pv) +% pv) ^ pv) | eq;
        var ph = mv | ~(xh | pv);
        var mh = pv & xh;
        if (ph & top != 0) score += 1;
        if (mh & top != 0) score -= 1;
        ph = (ph << 1) | 1;
        mh <<= 1;
        pv = mh | ~(xv | ph);
        mv = ph & xv;
    }
    return score;
}

pub const myers_ops = 15;
