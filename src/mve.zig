//! Minimal viable experiment: derive kernels, print them, verify 10^6 random
//! pairs against the scalar oracle, and time them against Myers' kernel.

const std = @import("std");
const bitdp = @import("bitdp");

fn printPlan(comptime s: bitdp.Scheme) void {
    const plan = comptime bitdp.derive(s);
    std.debug.print("\nscheme match={d} mismatch={d} gap={d} classes={d}: differences {any}, {d} ops/word, {d} carry chain(s), {s}\n", .{
        s.match, s.mismatch, s.gap, plan.cls.n, plan.vals[0..plan.k], plan.cost, plan.chains, @tagName(plan.method),
    });
    if (plan.len <= 40) for (plan.nodes[0..plan.len], 0..) |n, i| {
        std.debug.print("  r{d:<3} = {s:<5} {d} {d}\n", .{ i, @tagName(n.op), n.a, n.b });
    };
    std.debug.print("  out_v {any}  out_u {any}  lvl_pol {any}\n", .{ plan.out_v[0 .. plan.k - 1], plan.out_u[0 .. plan.k - 1], plan.lvl_pol[1..plan.k] });
}

fn verify(comptime s: bitdp.Scheme, pairs: usize) !void {
    const K = bitdp.Kernel(s);
    var prng = std.Random.DefaultPrng.init(0xb17d9);
    const rnd = prng.random();
    var buf: [65]i64 = undefined;
    var p: [64]u8 = undefined;
    var t: [256]u8 = undefined;
    for (0..pairs) |_| {
        const m = rnd.intRangeAtMost(usize, 1, 64);
        const n = rnd.intRangeAtMost(usize, 0, 256);
        for (p[0..m]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
        for (t[0..n]) |*x| x.* = s.alphabet[rnd.uintLessThan(usize, s.alphabet.len)];
        const want = bitdp.reference.scalar(s, p[0..m], t[0..n], &buf);
        const got = K.distance(p[0..m], t[0..n]);
        if (want != got) {
            std.debug.print("MISMATCH m={d} n={d} want={d} got={d}\n", .{ m, n, want, got });
            return error.Mismatch;
        }
    }
    std.debug.print("  verified {d} random pairs (m<=64, n<=256): bit-exact\n", .{pairs});
}

/// BLOSUM62 scores (NCBI, ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62), 20 standard residues.
const blosum62_order = "ARNDCQEGHILKMFPSTWYV";
const blosum62 = [20][20]i8{
    .{ 4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0 },
    .{ -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3 },
    .{ -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3 },
    .{ -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3 },
    .{ 0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },
    .{ -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2 },
    .{ -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2 },
    .{ 0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3 },
    .{ -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3 },
    .{ -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3 },
    .{ -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1 },
    .{ -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2 },
    .{ -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1 },
    .{ -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1 },
    .{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2 },
    .{ 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2 },
    .{ 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0 },
    .{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3 },
    .{ -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1 },
    .{ 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4 },
};

const blosum_index: [256]u8 = blk: {
    var t = [_]u8{0} ** 256;
    for (blosum62_order, 0..) |c, i| t[c] = i;
    break :blk t;
};

/// BLOSUM62 as a cost (negated score).
fn blosumCost(a: u8, b: u8) i32 {
    return -@as(i32, blosum62[blosum_index[a]][blosum_index[b]]);
}

/// DNA with transitions (A<->G, C<->T) cheaper than transversions.
fn tsTvCost(a: u8, b: u8) i32 {
    if (a == b) return 0;
    const purine = struct {
        fn f(x: u8) bool {
            return x == 'A' or x == 'G';
        }
    }.f;
    return if (purine(a) == purine(b)) 1 else 2;
}

const general = [_]bitdp.Scheme{
    .{ .sub = &tsTvCost, .gap = 2 },
    .{ .sub = &blosumCost, .alphabet = blosum62_order, .gap = 4 },
};

const bitpal = [_]bitdp.Scheme{
    .{ .match = 0, .mismatch = 1, .gap = 1 },
    .{ .match = -2, .mismatch = 3, .gap = 5 },
    .{ .match = -3, .mismatch = 4, .gap = 6 },
    .{ .match = -4, .mismatch = 5, .gap = 9 },
    .{ .match = -4, .mismatch = 7, .gap = 11 },
};

fn now(io: std.Io) i96 {
    return std.Io.Timestamp.now(io, .awake).nanoseconds;
}

/// Derived kernel vs the plain scalar DP (our unvectorized oracle), m = 64.
fn vsScalar(io: std.Io, gpa: std.mem.Allocator, comptime s: bitdp.Scheme, name: []const u8) !void {
    const n: usize = 2_000_000;
    const text = try gpa.alloc(u8, n);
    defer gpa.free(text);
    var prng = std.Random.DefaultPrng.init(11);
    for (text) |*x| x.* = s.alphabet[prng.random().uintLessThan(usize, s.alphabet.len)];
    const pat = text[5000..5064];
    var buf: [65]i64 = undefined;
    var t0 = now(io);
    const a = bitdp.reference.scalar(s, pat, text, &buf);
    const t_scalar = now(io) - t0;
    t0 = now(io);
    const b = bitdp.Kernel(s).distance(pat, text);
    const t_derived = now(io) - t0;
    if (a != b) return error.Mismatch;
    const cells: f64 = @floatFromInt(64 * n);
    std.debug.print("{s}: scalar {d:.2} GCUPS, derived {d:.2} GCUPS ({d} ops/word)\n", .{
        name, cells / @as(f64, @floatFromInt(t_scalar)), cells / @as(f64, @floatFromInt(t_derived)), bitdp.Kernel(s).ops,
    });
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    printPlan(.{});
    printPlan(.{ .mismatch = 2 });
    printPlan(.{ .mismatch = 1, .gap = 2 });
    printPlan(.{ .mismatch = 3, .gap = 2 });
    // BitPAl's published weight sets (score M, I, G), as costs (-M, -I, -G).
    inline for (bitpal) |w| printPlan(w);
    inline for (general) |w| printPlan(w);

    try verify(.{}, 1_000_000);
    try verify(.{ .mismatch = 2 }, 200_000);
    try verify(.{ .mismatch = 1, .gap = 2 }, 200_000);
    try verify(.{ .mismatch = 3, .gap = 2 }, 200_000);
    inline for (bitpal) |w| try verify(w, 100_000);
    inline for (general) |w| try verify(w, 100_000);

    // Throughput: m = 64 against a long random text.
    const gpa = init.gpa;
    const n: usize = 50_000_000;
    const text = try gpa.alloc(u8, n);
    defer gpa.free(text);
    var prng = std.Random.DefaultPrng.init(7);
    for (text) |*x| x.* = "ACGT"[prng.random().int(u2)];
    const pat = text[1000..1064];
    const K = bitdp.Kernel(.{});

    var t0 = now(io);
    const a = bitdp.reference.myers(pat, text);
    const t_myers = now(io) - t0;
    t0 = now(io);
    const b = K.distance(pat, text);
    const t_derived = now(io) - t0;
    std.mem.doNotOptimizeAway(a);
    std.mem.doNotOptimizeAway(b);
    if (a != b) return error.Mismatch;
    const cells: f64 = @floatFromInt(64 * n);
    try vsScalar(io, gpa, bitpal[1], "BitPAl (2,-3,-5)");
    try vsScalar(io, gpa, general[1], "BLOSUM62, gap 4");
    std.debug.print("\nthroughput m=64 n={d}: myers {d:.2} GCUPS, derived {d:.2} GCUPS (ops: myers {d}, derived {d})\n", .{
        n,
        cells / @as(f64, @floatFromInt(t_myers)),
        cells / @as(f64, @floatFromInt(t_derived)),
        bitdp.reference.myers_ops,
        K.ops,
    });
}
