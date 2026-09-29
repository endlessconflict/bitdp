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
//! a cascade of additions and bitwise logic. Two builders produce it: exact
//! truth-table synthesis for small difference sets, and a symbolic one that
//! decomposes the min-plus cell into thresholds and reads most planes off
//! unary sums computed by merging networks. Nothing is hand-derived.

const std = @import("std");

/// Costs to minimize. Row 0 and column 0 hold i * gap (global alignment).
pub const Scheme = struct {
    match: i32 = 0,
    mismatch: i32 = 1,
    gap: i32 = 1,
    /// Optional substitution cost for a (pattern, text) byte pair. When set,
    /// it replaces match/mismatch, and `alphabet` must list every byte that
    /// can occur so the compiler can collect the distinct costs.
    sub: ?*const fn (u8, u8) i32 = null,
    alphabet: []const u8 = "ACGT",
    /// `global`: the whole pattern against the whole text. `search`: the
    /// whole pattern against the best-matching substring of the text (row 0
    /// is all zeros and the result is the minimum over end positions).
    mode: Mode = .global,

    pub const Mode = enum { global, search };

    /// Horizontal difference in row 0: gap for global, 0 for search.
    pub fn topBoundary(s: Scheme) i32 {
        return if (s.mode == .search) 0 else s.gap;
    }

    pub fn cost(s: Scheme, a: u8, b: u8) i32 {
        if (s.sub) |f| return f(a, b);
        return if (a == b) s.match else s.mismatch;
    }
};

pub const max_classes = 24;
const max_inputs = max_vals + max_classes;

/// Distinct substitution costs, ascending. Class planes [cost >= c[j]] for
/// j >= 1 are the kernel's per-text-character inputs.
pub const Classes = struct { c: [max_classes]i32, n: u16 };

fn classes(s: Scheme) Classes {
    var r: Classes = .{ .c = undefined, .n = 0 };
    if (s.sub == null) {
        r.c[0] = @min(s.match, s.mismatch);
        r.c[1] = @max(s.match, s.mismatch);
        r.n = if (s.match == s.mismatch) 1 else 2;
        return r;
    }
    for (s.alphabet) |a| for (s.alphabet) |b| {
        const x = s.cost(a, b);
        if (std.mem.indexOfScalar(i32, r.c[0..r.n], x) == null) {
            if (r.n == max_classes) @compileError("bitdp: too many distinct substitution costs");
            r.c[r.n] = x;
            r.n += 1;
        }
    };
    std.mem.sort(i32, r.c[0..r.n], {}, std.sort.asc(i32));
    return r;
}

pub const max_vals = 32;
// ponytail: the truth-table path is exponential in variables; it only runs
// for k <= 5 distinct differences, where it sometimes beats the symbolic one.
const qm_max_vals = 5;
const max_vars = 2 * qm_max_vals - 1;
const TT = std.meta.Int(.unsigned, 1 << max_vars);
const max_nodes = 2048;
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
    table: [table_size]u16 = [_]u16{empty} ** table_size,

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
        // Hash-consing: identical nodes are shared (common subexpressions).
        var h: usize = (@as(usize, @intFromEnum(n.op)) *% 0x9e37 ^ @as(usize, n.a) *% 0x85eb ^ @as(usize, n.b) *% 0xc2b3) % table_size;
        while (p.table[h] != empty) : (h = (h + 1) % table_size) {
            const m = p.nodes[p.table[h]];
            if (m.op == n.op and m.a == n.a and m.b == n.b) return p.table[h];
        }
        if (p.len == max_nodes) @compileError("bitdp: program too large");
        p.nodes[p.len] = n;
        p.table[h] = p.len;
        p.len += 1;
        return p.len - 1;
    }
};

const table_size = 2 * max_nodes;
const empty = std.math.maxInt(u16);

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
fn cell(s: Scheme, dv_in: i32, sub: i32, dh_in: i32) Cell {
    const d = @min(sub, @min(dh_in + s.gap, dv_in + s.gap));
    return .{ .dv = d - dh_in, .dh = d - dv_in };
}

const Diffs = struct { vals: [max_vals]i32, k: u16 };

/// Every difference value the recurrence can produce (fixpoint from the
/// boundary). dv and dh obey the same recurrence and boundary, so one set.
fn diffs(s: Scheme, cls: Classes) Diffs {
    var set: [max_vals]i32 = undefined;
    set[0] = s.gap;
    var k: u16 = 1;
    if (s.topBoundary() != s.gap) {
        set[1] = s.topBoundary();
        k = 2;
    }
    var changed = true;
    while (changed) {
        changed = false;
        for (0..k) |i| for (0..k) |j| for (cls.c[0..cls.n]) |sub| {
            const c = cell(s, set[i], sub, set[j]);
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
    /// Kernel input a is stored complemented. Inputs 0..k-2 are the Delta-v
    /// threshold planes, then one plane per substitution class above the lowest.
    inp_pol: [max_inputs]bool,
    /// Level t is read complemented.
    lvl_pol: [max_vals]bool,
    /// Node computing next-column Delta-v threshold plane t (index t-1), in
    /// that input's stored polarity.
    out_v: [max_vals]u16,
    /// Node whose bit i holds level t of the state at row i-1 (index t-1).
    /// The kernel reads the bottom row's state at bit m.
    out_u: [max_vals]u16,
    /// Levels that became carry chains (one addition each).
    chains: u16,
    /// Word operations per column word, dead nodes excluded.
    cost: u32,
    /// Substitution cost classes; kernel inputs k-1.. are [cost >= c[j]], j >= 1.
    cls: Classes,
    /// Which builder produced the plan.
    method: enum { truth_table, direct, merge },
};

/// Input-variable view of one truth-table row.
const Row = struct { valid: bool, vi: u16, eq: bool };

fn decodeInput(a: u16, k: u16, in_pol: bool, cls_pol: bool) Row {
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
    return .{ .valid = valid, .vi = count, .eq = (a >> @intCast(k - 1) & 1 == 1) == cls_pol };
}

fn build(s: Scheme, d: Diffs, cls: Classes, in_pol: bool, cls_pol: bool) Plan {
    @setEvalBranchQuota(1 << 30);
    const k = d.k;
    const nl = k - 1; // levels t = 1..k-1

    // theta[t][vi][eq]: source of bit [next dh >= vals[t]].
    var theta: [max_vals][max_vals][2]Src = undefined;
    for (1..k) |t| for (0..k) |vi| for (0..2) |e| {
        var first: ?u16 = null;
        for (0..k) |si| {
            const g = cell(s, d.vals[vi], if (e == 1) cls.c[0] else cls.c[1], d.vals[si]).dh;
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
    plan.method = .truth_table;

    for (order[0..nl]) |t| {
        const top_bit = s.topBoundary() >= d.vals[t]; // row 0's dh
        var best: ?struct { p: Prog, prev: u16, u: u16, chain: bool, pol: bool } = null;
        for ([2]bool{ false, true }) |pol| {
            var q = p;
            var f_on: TT = 0;
            var p_on: TT = 0;
            var dc: TT = 0;
            for (0..@as(u16, 1) << @intCast(nv)) |ai| {
                const a: u16 = @intCast(ai);
                const row = decodeInput(a, k, in_pol, cls_pol);
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
        plan.out_u[t - 1] = best.?.prev;
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
            const row = decodeInput(a, k, in_pol, cls_pol);
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
            const dv = cell(s, d.vals[row.vi], if (row.eq) cls.c[0] else cls.c[1], d.vals[count]).dv;
            if ((dv >= d.vals[t]) != in_pol) on |= bit;
        }
        plan.out_v[t - 1] = synth(&p, nv, on, dc, vars[0..nv]);
    }

    return finish(p, plan, d, cls, cls_pol, in_pol);
}

/// Fill the common Plan fields and count live word operations.
fn finish(p: Prog, plan0: Plan, d: Diffs, cls: Classes, cls_pol: bool, in_pol: bool) Plan {
    var plan = plan0;
    plan.vals = d.vals;
    plan.k = d.k;
    plan.nodes = p.nodes;
    plan.len = p.len;
    plan.cls = cls;
    for (0..d.k - 1) |a| plan.inp_pol[a] = in_pol;
    for (0..cls.n -| 1) |j| plan.inp_pol[d.k - 1 + j] = cls_pol;
    plan.cost = liveCost(plan);
    return plan;
}

const Planes = enum { v, u };
const max_unary = 64;
const Unary = struct { bits: [2 * max_unary]u16 = undefined, len: usize };

/// Word operations of the nodes the outputs actually use.
fn liveCost(plan: Plan) u32 {
    const live = liveNodes(plan);
    var cost: u32 = 0;
    for (plan.nodes[0..plan.len], 0..) |n, j| {
        if (live[j]) cost += opCost(n.op);
    }
    return cost;
}

fn liveNodes(plan: Plan) [max_nodes]bool {
    const nl = plan.k - 1;
    var live = [_]bool{false} ** max_nodes;
    for (plan.out_v[0..nl]) |o| live[o] = true;
    for (plan.out_u[0..nl]) |o| live[o] = true;
    var i: usize = plan.len;
    while (i > 0) {
        i -= 1;
        if (!live[i]) continue;
        const n = plan.nodes[i];
        switch (n.op) {
            .input, .zero, .ones => {},
            .not, .shl0, .shl1 => live[n.a] = true,
            else => {
                live[n.a] = true;
                live[n.b] = true;
            },
        }
    }
    return live;
}

/// Threshold planes and the operations the symbolic builder composes.
const Sym = struct {
    p: *Prog,
    s: Scheme,
    vals: []const i32,
    in_pol: bool,
    lvl_pol: bool,
    xin: [max_vals]u16 = undefined,
    uprev: [max_vals]u16 = undefined,
    zero: u16,
    ones: u16,
    cls: Classes,
    cls_pol: bool,
    /// cpl[j] = [cost >= cls.c[j]] for j >= 1, stored complemented if cls_pol.
    cpl: [max_classes]u16 = undefined,

    fn op(y: *Sym, o: Op, a: u16, b: u16) u16 {
        return y.p.emit(.{ .op = o, .a = a, .b = b });
    }

    /// [x >= w] (positive) or its complement, over the Delta-v input planes
    /// (`v`) or the previous-row state planes (`u`).
    fn ge(y: *Sym, which: Planes, w: i32, positive: bool) u16 {
        const k = y.vals.len;
        if (w <= y.vals[0]) return if (positive) y.ones else y.zero;
        if (w > y.vals[k - 1]) return if (positive) y.zero else y.ones;
        var t: usize = 1;
        while (y.vals[t] < w) t += 1;
        const stored = if (which == .v) y.xin[t] else y.uprev[t];
        const pol = if (which == .v) y.in_pol else y.lvl_pol;
        return if (positive != pol) stored else y.op(.not, stored, 0);
    }

    /// [c_e - x >= t] = not [x >= c_e - t + 1], with c_e the match or
    /// mismatch cost. match <= mismatch makes the match side a subset.
    /// Batcher's odd-even merge of two descending 0/1 sequences of the same
    /// power-of-two length. A comparator on bits is (x or y, x and y).
    fn merge(y: *Sym, a: []const u16, b: []const u16, out: []u16) void {
        const n = a.len;
        if (n == 1) {
            out[0] = y.op(.@"or", a[0], b[0]);
            out[1] = y.op(.@"and", a[0], b[0]);
            return;
        }
        var ae: [max_unary]u16 = undefined;
        var ao: [max_unary]u16 = undefined;
        var be: [max_unary]u16 = undefined;
        var bo: [max_unary]u16 = undefined;
        for (0..n / 2) |i| {
            ae[i] = a[2 * i];
            ao[i] = a[2 * i + 1];
            be[i] = b[2 * i];
            bo[i] = b[2 * i + 1];
        }
        var v: [2 * max_unary]u16 = undefined;
        var w: [2 * max_unary]u16 = undefined;
        y.merge(ae[0 .. n / 2], be[0 .. n / 2], v[0..n]);
        y.merge(ao[0 .. n / 2], bo[0 .. n / 2], w[0..n]);
        out[0] = v[0];
        for (0..n - 1) |i| {
            out[2 * i + 1] = y.op(.@"or", w[i], v[i + 1]);
            out[2 * i + 2] = y.op(.@"and", w[i], v[i + 1]);
        }
        out[2 * n - 1] = w[n - 1];
    }

    /// Unary addition: given thermometer codes of two counts, the merged
    /// sequence is the thermometer code of their sum, bit t-1 = [sum >= t].
    fn unarySum(y: *Sym, a: []const u16, b: []const u16) Unary {
        var n: usize = 1;
        while (n < a.len or n < b.len) n *= 2;
        var pa = [_]u16{y.zero} ** max_unary;
        var pb = [_]u16{y.zero} ** max_unary;
        @memcpy(pa[0..a.len], a);
        @memcpy(pb[0..b.len], b);
        var r: Unary = .{ .len = a.len + b.len };
        y.merge(pa[0..n], pb[0..n], r.bits[0 .. 2 * n]);
        return r;
    }

    fn unaryGe(y: *Sym, z: Unary, t: i32) u16 {
        if (t <= 0) return y.ones;
        if (t > z.len) return y.zero;
        return z.bits[@intCast(t - 1)];
    }

    /// [c_e >= x] and plane: all of it, only mismatch rows, or nothing.
    fn costMask(y: *Sym, x: i32, plane: u16) u16 {
        return y.op(.@"and", y.costAtLeast(x), plane);
    }

    /// [c_e >= x] as a class plane or a constant.
    fn costAtLeast(y: *Sym, x: i32) u16 {
        const c = y.cls.c[0..y.cls.n];
        if (x <= c[0]) return y.ones;
        if (x > c[c.len - 1]) return y.zero;
        var j: usize = 1;
        while (c[j] < x) j += 1;
        return y.maybeNot(y.cpl[j], y.cls_pol);
    }

    fn maybeNot(y: *Sym, x: u16, complement: bool) u16 {
        return if (complement) y.op(.not, x, 0) else x;
    }

    fn costGe(y: *Sym, which: Planes, t: i32) u16 {
        // [c_e - x >= t] = OR over classes j of [c_e >= c_j] and [x <= c_j - t].
        var r = y.zero;
        for (y.cls.c[0..y.cls.n]) |c| r = y.op(.@"or", r, y.op(.@"and", y.costAtLeast(c), y.ge(which, c - t + 1, false)));
        return r;
    }
};

/// Symbolic derivation by threshold decomposition of the min-plus cell
/// d = min(c_e, s + gap, v + gap), using
///   [min(a, b) >= t] = [a >= t] and [b >= t],
///   [x - y >= q]     = OR over w of [x >= w] and not [y >= w - q + 1].
/// Every level and output plane becomes a short OR of ANDs over threshold
/// planes. Size grows like k^2 rather than 2^k, so wide score ranges work.
fn buildSym(s: Scheme, d: Diffs, cls: Classes, cls_pol: bool, in_pol: bool, lvl_pol: bool, bulk: bool) Plan {
    @setEvalBranchQuota(1 << 30);
    const k = d.k;
    const vals = d.vals[0..k];
    var p = Prog{};
    var plan: Plan = undefined;
    plan.chains = 0;
    plan.method = if (bulk) .merge else .direct;
    const zero = p.emit(.{ .op = .zero });
    const ones = p.emit(.{ .op = .ones });
    var y = Sym{ .p = &p, .s = s, .vals = vals, .in_pol = in_pol, .lvl_pol = lvl_pol, .zero = zero, .ones = ones, .cls = cls, .cls_pol = cls_pol };
    for (1..k) |t| y.xin[t] = p.emit(.{ .op = .input, .a = @intCast(t - 1) });
    for (1..cls.n) |j| y.cpl[j] = p.emit(.{ .op = .input, .a = @intCast(k - 1 + j - 1) });
    const cmin = cls.c[0];
    const cmax = cls.c[cls.n - 1];
    const vmin = vals[0];
    const vmax = vals[k - 1];
    // Every level's terms stop depending on v once the source level reaches
    // mismatch - gap (see below), independently of the level. Levels above
    // that cap need no carry chain, and in bulk mode all of them are read
    // off one unary sum: [min(s, cap) - v >= q] for every q at once.
    const cap = @max(cmax - s.gap, vmin);
    var zlev: ?Unary = null;
    var zlev_lo: i32 = 0;

    // Level t: u_t = [gap >= t] and [c_e - v >= t] and [s - v >= t - gap].
    // Terms of the last factor with source level r > t are absorbed, and the
    // r = t term has condition [v <= gap], which always holds, so
    //   u_t = C and (Q or u_t(i-1)),  C = [gap >= t] and [c_e - v >= t],
    // with Q the terms of levels r < t. That is a carry chain with generate
    // F = C and Q, propagating wherever C holds.
    for (1..k) |t| {
        const tau = vals[t];
        const cin = s.topBoundary() >= tau; // row 0's dh
        const c = if (s.gap >= tau) y.costGe(.v, tau) else zero;
        const q = tau - s.gap;
        // C implies v <= mismatch - tau. Once a term's own bound on v reaches
        // that, its v-condition is implied by C: the term becomes the bare
        // plane u_r, which contains every later term and the self term u_t.
        // Then the level needs no carry chain at all.
        const bound = cmax - tau;
        var covered = vals[0] - q >= bound;
        var qn = if (covered) ones else y.ge(.v, vals[0] - q + 1, false);
        if (!covered) for (1..t) |r| {
            if (vals[r] - q >= bound) {
                qn = y.op(.@"or", qn, y.ge(.u, vals[r], true));
                covered = true;
                break;
            }
            qn = y.op(.@"or", qn, y.op(.@"and", y.ge(.u, vals[r], true), y.ge(.v, vals[r] - q + 1, false)));
        };
        var f = y.op(.@"and", c, qn);
        if (bulk and tau > cap) {
            // dh_out = min(w' - v, gap) with w' = min(c_e, min(s, cap) + gap),
            // so the level is one bit of the unary sum of w' and -v.
            if (zlev == null) {
                const lo: i32 = @min(cmin, vmin + s.gap);
                const hi: i32 = @min(cmax, cap + s.gap);
                var a: [max_unary]u16 = undefined;
                var b: [max_unary]u16 = undefined;
                const na: usize = @intCast(hi - lo);
                const nb: usize = @intCast(vmax - vmin);
                for (0..na) |j| {
                    const x = lo + @as(i32, @intCast(j)) + 1;
                    const sx = if (x - s.gap <= cap) y.ge(.u, x - s.gap, true) else zero;
                    a[j] = y.costMask(x, sx);
                }
                for (0..nb) |j| b[j] = y.ge(.v, vmax - @as(i32, @intCast(j)), false);
                zlev = y.unarySum(a[0..na], b[0..nb]);
                zlev_lo = lo;
            }
            covered = true;
            f = y.unaryGe(zlev.?, tau - zlev_lo + vmax);
        }
        var prev: u16 = undefined;
        if (covered or qn == ones or c == zero) {
            const c_or_f = if (covered) f else c;
            const u = if (lvl_pol) y.op(.not, c_or_f, 0) else c_or_f;
            prev = y.op(if (cin != lvl_pol) .shl1 else .shl0, u, 0);
        } else {
            // a = C and g = F (complement: a = ~F, g = ~C). The carries of
            // a + g + cin are u(i-1), and sum ^ a ^ g = sum ^ (C ^ F).
            const a = if (lvl_pol) y.op(.not, f, 0) else c;
            const g = if (lvl_pol) y.op(.not, c, 0) else f;
            const sum = y.op(if (cin != lvl_pol) .add1 else .add, a, g);
            prev = y.op(.xor, sum, y.op(.xor, c, f));
            plan.chains += 1;
        }
        plan.lvl_pol[t] = lvl_pol;
        plan.out_u[t - 1] = prev;
        y.uprev[t] = prev;
    }

    // Output: [dv_out >= t] = [gap >= t] and [c_e - s >= t] and [v - s >= t - gap].
    var zout: ?Unary = null;
    var zout_lo: i32 = 0;
    for (1..k) |t| {
        const tau = vals[t];
        var out = zero;
        if (s.gap >= tau) {
            const q = tau - s.gap;
            // [c_e - s >= t] implies s <= mismatch - t, so once w - q reaches
            // that bound the s-condition is implied: the term is [v >= w]
            // alone, and it contains every later term.
            const bound = cmax - tau;
            var r = zero;
            if (bulk) {
                // dv_out = min(w - s, gap) with w = min(c_e, v + gap), which
                // depends on inputs only: one unary sum of w and -s gives
                // every output plane, with no per-plane mask.
                if (zout == null) {
                    const lo: i32 = @min(cmin, vmin + s.gap);
                    const hi: i32 = @min(cmax, vmax + s.gap);
                    var a: [max_unary]u16 = undefined;
                    var b: [max_unary]u16 = undefined;
                    const na: usize = @intCast(hi - lo);
                    const nb: usize = @intCast(vmax - vmin);
                    for (0..na) |j| {
                        const x = lo + @as(i32, @intCast(j)) + 1;
                        a[j] = y.costMask(x, y.ge(.v, x - s.gap, true));
                    }
                    for (0..nb) |j| b[j] = y.ge(.u, vmax - @as(i32, @intCast(j)), false);
                    zout = y.unarySum(a[0..na], b[0..nb]);
                    zout_lo = lo;
                }
                plan.out_v[t - 1] = y.maybeNot(y.unaryGe(zout.?, tau - zout_lo + vmax), in_pol);
                continue;
            } else for (vals) |w| {
                if (w - q >= bound) {
                    r = y.op(.@"or", r, y.ge(.v, w, true));
                    break;
                }
                r = y.op(.@"or", r, y.op(.@"and", y.ge(.v, w, true), y.ge(.u, w - q + 1, false)));
            }
            out = y.op(.@"and", y.costGe(.u, tau), r);
        }
        plan.out_v[t - 1] = if (in_pol) y.op(.not, out, 0) else out;
    }
    return finish(p, plan, d, cls, cls_pol, in_pol);
}

/// Which builders and passes `deriveWith` may use (for ablations).
pub const Options = struct {
    /// Exhaustive truth-table synthesis, tried when k <= 5 and two cost classes.
    truth_table: bool = true,
    /// Symbolic builder, every threshold as its own OR of ANDs.
    direct: bool = true,
    /// Symbolic builder with unary sums from merging networks.
    merge: bool = true,
};

/// Derive the cheapest plan with every builder and pass enabled.
pub fn derive(comptime s: Scheme) Plan {
    return deriveWith(s, .{});
}

/// Derive the cheapest plan from the enabled builders.
pub fn deriveWith(comptime s: Scheme, comptime opt: Options) Plan {
    comptime {
        @setEvalBranchQuota(1 << 30);
        if (s.gap <= 0) @compileError("bitdp: gap cost must be positive");
        const cls = classes(s);
        const d = diffs(s, cls);
        var cand: [24]Plan = undefined;
        var nc: usize = 0;
        // ponytail: one plane polarity per builder. Building every polarity
        // combination multiplied compile memory for no op-count gain, and a
        // single-flip De Morgan search never removed a NOT (NOTs sit on values
        // needed in both polarities; removing one means flipping a whole
        // subnetwork). A coordinated search is the upgrade path.
        for ([2]bool{ false, true }) |bulk| {
            if (if (bulk) !opt.merge else !opt.direct) continue;
            cand[nc] = buildSym(s, d, cls, false, false, false, bulk);
            nc += 1;
        }
        if (opt.truth_table and d.k <= qm_max_vals and cls.n == 2) for (0..4) |pol| {
            cand[nc] = build(s, d, cls, pol & 1 != 0, pol & 2 != 0);
            nc += 1;
        };
        if (nc == 0) @compileError("bitdp: no builder enabled for this scheme");
        var best = cand[0];
        for (cand[1..nc]) |c| {
            if (c.cost < best.cost) best = c;
        }
        return best;
    }
}
