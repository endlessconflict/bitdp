//! Compile-time derivation of a bit-parallel column update for a DP scheme.
//!
//! The DP is read one text column at a time. Pattern rows live in the bits
//! of a machine word. Inside a column, the value passed from row i-1 to row i
//! is the horizontal difference dh. With bounded differences this makes the
//! column a finite-state machine whose state is dh. Its step functions are
//! monotone in dh, so every threshold bit [dh >= t] of the next state is a
//! constant or a copy of one threshold bit of the previous state. That splits
//! the machine into one "level" per threshold:
//!
//!   * a level that copies itself is a carry chain (generate / propagate /
//!     kill), and a carry chain over a word is exactly one integer addition;
//!   * a level that copies another level is a shift plus bitwise logic.
//!
//! If the copy graph between distinct levels is acyclic, the whole column is
//! a cascade of additions and bitwise logic. Every boolean function in that
//! cascade is synthesized here from its truth table. Nothing is hand-derived.

const std = @import("std");

/// Costs to minimize. Row 0 and column 0 hold i * gap (global alignment).
pub const Scheme = struct {
    match: i32 = 0,
    mismatch: i32 = 1,
    gap: i32 = 1,
};

// ponytail: fixed caps keep comptime cheap; 2k-1 <= 9 variables means up to
// 5 distinct differences. Wider score ranges need a smarter synthesizer.
pub const max_vals = 5;
const max_vars = 2 * max_vals - 1;
const TT = std.meta.Int(.unsigned, 1 << max_vars);
const max_nodes = 256;
const max_primes = 512;

pub const Op = enum(u8) { input, zero, ones, not, @"and", @"or", xor, add, add1, shl0, shl1 };

/// `input` reads variable `a`. `add1` is a + b + 1. `shl0`/`shl1` shift left
/// by one and fill bit 0 with 0/1.
pub const Node = struct { op: Op, a: u16 = 0, b: u16 = 0 };

/// Word operations a node costs. a + b + 1 and (x << 1) | 1 count as two.
pub fn opCost(op: Op) u32 {
    return switch (op) {
        .input, .zero, .ones => 0,
        .add1, .shl1 => 2,
        else => 1,
    };
}

const Prog = struct {
    nodes: [max_nodes]Node = undefined,
    len: u16 = 0,

    fn emit(p: *Prog, n0: Node) u16 {
        var n = n0;
        const binary = switch (n.op) {
            .@"and", .@"or", .xor, .add, .add1 => true,
            else => false,
        };
        if (binary and n.a > n.b) std.mem.swap(u16, &n.a, &n.b);
        if (n.op == .not) switch (p.nodes[n.a].op) {
            .not => return p.nodes[n.a].a,
            .zero => return p.emit(.{ .op = .ones }),
            .ones => return p.emit(.{ .op = .zero }),
            else => {},
        };
        if (n.op == .@"and" or n.op == .@"or") {
            if (n.a == n.b) return n.a;
            const absorbing: Op = if (n.op == .@"and") .zero else .ones;
            const neutral: Op = if (n.op == .@"and") .ones else .zero;
            for ([2]u16{ n.a, n.b }, [2]u16{ n.b, n.a }) |x, y| {
                if (p.nodes[x].op == absorbing) return x;
                if (p.nodes[x].op == neutral) return y;
            }
        }
        for (p.nodes[0..p.len], 0..) |m, i| {
            if (m.op == n.op and m.a == n.a and m.b == n.b) return @intCast(i);
        }
        p.nodes[p.len] = n;
        p.len += 1;
        return p.len - 1;
    }
};

// ---------------------------------------------------------------- synthesis

fn fullMask(nv: u16) TT {
    const size: u16 = @as(u16, 1) << @intCast(nv);
    return if (size == @bitSizeOf(TT)) ~@as(TT, 0) else (@as(TT, 1) << @intCast(size)) - 1;
}

/// Minterms covered by the cube {x : x & care == val}.
fn cubeCover(nv: u16, care: u16, val: u16) TT {
    const free = ~care & ((@as(u16, 1) << @intCast(nv)) - 1);
    var cover: TT = 0;
    var s: u16 = 0;
    while (true) {
        cover |= @as(TT, 1) << @intCast(val | s);
        if (s == free) break;
        s = (s -% free) & free;
    }
    return cover;
}

/// Sum of products: prime implicants + greedy cover (Quine-McCluskey style).
fn sop(p: *Prog, nv: u16, on: TT, dc: TT, vars: []const u16) u16 {
    const allowed = on | dc;
    var primes: [max_primes][2]u16 = undefined;
    var np: usize = 0;
    var care: u16 = 0;
    while (care < (@as(u16, 1) << @intCast(nv))) : (care += 1) {
        var val: u16 = care;
        while (true) {
            if (cubeCover(nv, care, val) & ~allowed == 0) {
                var prime = true;
                for (0..nv) |i| {
                    const bit = @as(u16, 1) << @intCast(i);
                    if (care & bit != 0 and cubeCover(nv, care & ~bit, val & ~bit) & ~allowed == 0) {
                        prime = false;
                        break;
                    }
                }
                if (prime) {
                    if (np == max_primes) @compileError("bitdp: too many prime implicants");
                    primes[np] = .{ care, val };
                    np += 1;
                }
            }
            if (val == 0) break;
            val = (val - 1) & care;
        }
    }

    var uncovered = on;
    var result: ?u16 = null;
    while (uncovered != 0) {
        var best: usize = 0;
        var best_gain: u32 = 0;
        for (primes[0..np], 0..) |c, i| {
            const gain: u32 = @popCount(cubeCover(nv, c[0], c[1]) & uncovered);
            if (gain > best_gain or (gain == best_gain and gain > 0 and
                @popCount(c[0]) < @popCount(primes[best][0])))
            {
                best = i;
                best_gain = gain;
            }
        }
        const c = primes[best];
        uncovered &= ~cubeCover(nv, c[0], c[1]);
        var term: ?u16 = null;
        for (0..nv) |i| {
            if (c[0] >> @intCast(i) & 1 == 0) continue;
            const lit = if (c[1] >> @intCast(i) & 1 == 1) vars[i] else p.emit(.{ .op = .not, .a = vars[i] });
            term = if (term) |t| p.emit(.{ .op = .@"and", .a = t, .b = lit }) else lit;
        }
        const t = term orelse p.emit(.{ .op = .ones });
        result = if (result) |r| p.emit(.{ .op = .@"or", .a = r, .b = t }) else t;
    }
    return result.?;
}

/// Cheapest of f and ~f (as a NOT of the complement's SOP).
fn synth(p: *Prog, nv: u16, on: TT, dc: TT, vars: []const u16) u16 {
    const off = fullMask(nv) & ~(on | dc);
    if (on == 0) return p.emit(.{ .op = .zero });
    if (off == 0) return p.emit(.{ .op = .ones });
    var a = p.*;
    const na = sop(&a, nv, on, dc, vars);
    var b = p.*;
    const nb = b.emit(.{ .op = .not, .a = sop(&b, nv, off, dc, vars) });
    if (b.len < a.len) {
        p.* = b;
        return nb;
    }
    p.* = a;
    return na;
}

// ---------------------------------------------------------------- the DP

const Cell = struct { dv: i32, dh: i32 };

/// One DP cell with the diagonal value taken as 0: left = dv_in, up = dh_in.
fn cell(s: Scheme, dv_in: i32, eq: bool, dh_in: i32) Cell {
    const sub = if (eq) s.match else s.mismatch;
    const d = @min(sub, @min(dh_in + s.gap, dv_in + s.gap));
    return .{ .dv = d - dh_in, .dh = d - dv_in };
}

const Diffs = struct { vals: [max_vals]i32, k: u16 };

/// Every difference value the recurrence can produce (fixpoint from the
/// boundary). dv and dh obey the same recurrence and boundary, so one set.
fn diffs(s: Scheme) Diffs {
    var set: [max_vals]i32 = undefined;
    set[0] = s.gap;
    var k: u16 = 1;
    var changed = true;
    while (changed) {
        changed = false;
        for (0..k) |i| for (0..k) |j| for ([2]bool{ false, true }) |eq| {
            const c = cell(s, set[i], eq, set[j]);
            for ([2]i32{ c.dv, c.dh }) |x| {
                if (std.mem.indexOfScalar(i32, set[0..k], x) == null) {
                    if (k == max_vals) @compileError("bitdp: score scheme has too many distinct differences");
                    set[k] = x;
                    k += 1;
                    changed = true;
                }
            }
        };
    }
    std.mem.sort(i32, set[0..k], {}, std.sort.asc(i32));
    return .{ .vals = set, .k = k };
}

/// Where threshold bit t of the next state comes from, for one input symbol.
const Src = union(enum) { zero, one, level: u16 };

pub const Plan = struct {
    vals: [max_vals]i32,
    k: u16,
    nodes: [max_nodes]Node,
    len: u16,
    /// Delta-v planes are stored complemented.
    in_pol: bool,
    /// Level t is carried complemented.
    lvl_pol: [max_vals]bool,
    /// Node computing next-column Delta-v threshold plane t (index t-1).
    out_v: [max_vals]u16,
    /// Node computing level t of the state at each row (index t-1).
    out_u: [max_vals]u16,
    /// Levels that became carry chains (one addition each).
    chains: u16,
    /// Word operations per column word, dead nodes excluded.
    cost: u32,
};

/// Input-variable view of one truth-table row.
const Row = struct { valid: bool, vi: u16, eq: bool };

fn decodeInput(a: u16, k: u16, in_pol: bool) Row {
    var count: u16 = 0;
    var seen_zero = false;
    var valid = true;
    for (1..k) |t| {
        const bit = (a >> @intCast(t - 1) & 1 == 1) != in_pol;
        if (bit) {
            if (seen_zero) valid = false;
            count += 1;
        } else seen_zero = true;
    }
    return .{ .valid = valid, .vi = count, .eq = a >> @intCast(k - 1) & 1 == 1 };
}

fn build(s: Scheme, d: Diffs, in_pol: bool) Plan {
    @setEvalBranchQuota(1 << 30);
    const k = d.k;
    const nl = k - 1; // levels t = 1..k-1

    // theta[t][vi][eq]: source of bit [next dh >= vals[t]].
    var theta: [max_vals][max_vals][2]Src = undefined;
    for (1..k) |t| for (0..k) |vi| for (0..2) |e| {
        var first: ?u16 = null;
        for (0..k) |si| {
            const g = cell(s, d.vals[vi], e == 1, d.vals[si]).dh;
            const bit = g >= d.vals[t];
            if (bit and first == null) first = @intCast(si);
            if (!bit and first != null) @compileError("bitdp: step function is not monotone");
        }
        theta[t][vi][e] = if (first) |f| (if (f == 0) .one else .{ .level = f }) else .zero;
    };

    // Topological order of levels over cross-level copies.
    var order: [max_vals]u16 = undefined;
    var done = [_]bool{false} ** max_vals;
    for (0..nl) |pos| {
        var pick: ?u16 = null;
        for (1..k) |t| {
            if (done[t]) continue;
            var ready = true;
            for (0..k) |vi| for (0..2) |e| switch (theta[t][vi][e]) {
                .level => |r| if (r != t and !done[r]) {
                    ready = false;
                },
                else => {},
            };
            if (ready) {
                pick = @intCast(t);
                break;
            }
        }
        // ponytail: cyclic copies need a general transition-monoid scan;
        // not implemented until a real scheme needs it.
        const t = pick orelse @compileError("bitdp: cyclic level dependencies not supported yet");
        order[pos] = t;
        done[t] = true;
    }

    var p = Prog{};
    var vars: [max_vars]u16 = undefined;
    for (0..k) |i| vars[i] = p.emit(.{ .op = .input, .a = @intCast(i) });
    var nv: u16 = k; // inputs: k-1 Delta-v planes + eq
    var var_of = [_]u16{0} ** max_vals; // level -> variable index of prev(u)
    var plan: Plan = undefined;
    plan.chains = 0;

    for (order[0..nl]) |t| {
        const top_bit = s.gap >= d.vals[t]; // row 0 holds dh = gap
        var best: ?struct { p: Prog, prev: u16, u: u16, chain: bool, pol: bool } = null;
        for ([2]bool{ false, true }) |pol| {
            var q = p;
            var f_on: TT = 0;
            var p_on: TT = 0;
            var dc: TT = 0;
            for (0..@as(u16, 1) << @intCast(nv)) |ai| {
                const a: u16 = @intCast(ai);
                const row = decodeInput(a, k, in_pol);
                const bit = @as(TT, 1) << @intCast(a);
                if (!row.valid) {
                    dc |= bit;
                    continue;
                }
                const f: bool = switch (theta[t][row.vi][@intFromBool(row.eq)]) {
                    .one => !pol,
                    .zero => pol,
                    .level => |r| if (r == t) blk: {
                        p_on |= bit;
                        break :blk false;
                    } else ((a >> @intCast(var_of[r]) & 1 == 1) != plan.lvl_pol[r]) != pol,
                };
                if (f) f_on |= bit;
            }
            const cin = top_bit != pol;
            const fn_ = synth(&q, nv, f_on, dc, vars[0..nv]);
            var prev: u16 = undefined;
            var u: u16 = undefined;
            if (p_on != 0) {
                // u_i = F_i | P_i & u_{i-1}. With a = F | P, the carries of
                // a + F + cin are exactly u_{i-1}, and sum ^ a ^ F = sum ^ P.
                const pn = synth(&q, nv, p_on, dc, vars[0..nv]);
                const an = synth(&q, nv, f_on | p_on, dc, vars[0..nv]);
                const sum = q.emit(.{ .op = if (cin) .add1 else .add, .a = an, .b = fn_ });
                prev = q.emit(.{ .op = .xor, .a = sum, .b = pn });
                u = q.emit(.{ .op = .@"or", .a = fn_, .b = q.emit(.{ .op = .@"and", .a = pn, .b = prev }) });
            } else {
                u = fn_;
                prev = q.emit(.{ .op = if (cin) .shl1 else .shl0, .a = fn_ });
            }
            if (best == null or q.len < best.?.p.len)
                best = .{ .p = q, .prev = prev, .u = u, .chain = p_on != 0, .pol = pol };
        }
        p = best.?.p;
        plan.lvl_pol[t] = best.?.pol;
        plan.out_u[t - 1] = best.?.u;
        if (best.?.chain) plan.chains += 1;
        vars[nv] = best.?.prev;
        var_of[t] = nv;
        nv += 1;
    }

    // Next-column Delta-v planes from inputs and every prev(u) level.
    for (1..k) |t| {
        var on: TT = 0;
        var dc: TT = 0;
        for (0..@as(u16, 1) << @intCast(nv)) |ai| {
            const a: u16 = @intCast(ai);
            const bit = @as(TT, 1) << @intCast(a);
            const row = decodeInput(a, k, in_pol);
            var count: u16 = 0;
            var seen_zero = false;
            var valid = row.valid;
            for (1..k) |l| {
                const ub = (a >> @intCast(var_of[l]) & 1 == 1) != plan.lvl_pol[l];
                if (ub) {
                    if (seen_zero) valid = false;
                    count += 1;
                } else seen_zero = true;
            }
            if (!valid) {
                dc |= bit;
                continue;
            }
            const dv = cell(s, d.vals[row.vi], row.eq, d.vals[count]).dv;
            if ((dv >= d.vals[t]) != in_pol) on |= bit;
        }
        plan.out_v[t - 1] = synth(&p, nv, on, dc, vars[0..nv]);
    }

    // Cost of live nodes only.
    var live = [_]bool{false} ** max_nodes;
    for (plan.out_v[0..nl]) |o| live[o] = true;
    for (plan.out_u[0..nl]) |o| live[o] = true;
    var i: usize = p.len;
    while (i > 0) {
        i -= 1;
        if (!live[i]) continue;
        const n = p.nodes[i];
        switch (n.op) {
            .input, .zero, .ones => {},
            .not, .shl0, .shl1 => live[n.a] = true,
            else => {
                live[n.a] = true;
                live[n.b] = true;
            },
        }
    }
    var cost: u32 = 0;
    for (p.nodes[0..p.len], 0..) |n, j| {
        if (live[j]) cost += opCost(n.op);
    }

    plan.vals = d.vals;
    plan.k = k;
    plan.nodes = p.nodes;
    plan.len = p.len;
    plan.in_pol = in_pol;
    plan.cost = cost;
    return plan;
}

/// Derive the cheapest plan over both Delta-v plane polarities.
pub fn derive(comptime s: Scheme) Plan {
    comptime {
        @setEvalBranchQuota(1 << 30);
        if (s.gap <= 0) @compileError("bitdp: gap cost must be positive");
        const d = diffs(s);
        const a = build(s, d, false);
        const b = build(s, d, true);
        return if (b.cost < a.cost) b else a;
    }
}
