//! Demo on a real genome: CRISPR-style 20-nt guides scanned along a whole
//! genome in search mode with transition/transversion costs, reporting every
//! end position with cost <= K. Hits are checked against the scalar DP.
//!
//! usage: bitdp-scan GENOME.fa GUIDES K [verify] [packed]
//!
//! With "packed", three guides share each lane word (Kernel.Packed).

const std = @import("std");
const bitdp = @import("bitdp");

const scheme: bitdp.Scheme = .{ .sub = &bitdp.schemes.tsTvCost, .gap = 2, .mode = .search };
const Kn = bitdp.Kernel(scheme);
const guide_len = 20;

fn now(io: std.Io) i96 {
    return std.Io.Timestamp.now(io, .awake).nanoseconds;
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.arena.allocator();
    const args = try init.minimal.args.toSlice(gpa);
    if (args.len < 4) return error.Usage;
    const raw = try std.Io.Dir.cwd().readFileAlloc(io, args[1], gpa, .unlimited);
    const nguides = try std.fmt.parseInt(usize, args[2], 10);
    const max_cost = try std.fmt.parseInt(i64, args[3], 10);
    var verify = false;
    var pack = false;
    for (args[4..]) |a| {
        if (std.mem.eql(u8, a, "verify")) verify = true;
        if (std.mem.eql(u8, a, "packed")) pack = true;
    }

    // FASTA: drop header lines and newlines, upper-case.
    var genome: std.ArrayList(u8) = .empty;
    var lines = std.mem.tokenizeScalar(u8, raw, '\n');
    while (lines.next()) |l| {
        if (l.len == 0 or l[0] == '>') continue;
        for (std.mem.trimEnd(u8, l, "\r")) |c| try genome.append(gpa, std.ascii.toUpper(c));
    }
    const g = genome.items;

    // Guides: genome 20-mers with two random substitutions, so each has a
    // near-exact site plus whatever the genome offers by chance.
    var prng = std.Random.DefaultPrng.init(2026);
    const rnd = prng.random();
    const guides = try gpa.alloc([guide_len]u8, nguides);
    for (guides) |*gd| {
        const at = rnd.uintLessThan(usize, g.len - guide_len);
        @memcpy(gd, g[at..][0..guide_len]);
        for (0..2) |_| gd[rnd.uintLessThan(usize, guide_len)] = "ACGT"[rnd.int(u2)];
    }

    var hits: std.ArrayList(Kn.Hit) = .empty;
    var total_hits: usize = 0;
    var by_guide = try gpa.alloc(usize, nguides);
    @memset(by_guide, 0);
    const t0 = now(io);
    var start: usize = 0;
    const slices = try gpa.alloc([]const u8, nguides);
    for (slices, guides) |*sl, *gd| sl.* = gd;
    while (pack and start < nguides) {
        const pk = Kn.Packed.init(slices[start..]);
        hits.clearRetainingCapacity();
        try pk.scan(gpa, g, max_cost, &hits);
        for (hits.items) |h| by_guide[start + h.lane] += 1;
        start += pk.count;
    }
    while (start < nguides) : (start += Kn.lanes) {
        const count = @min(Kn.lanes, nguides - start);
        var ps: [Kn.lanes][]const u8 = undefined;
        for (0..count) |l| ps[l] = &guides[start + l];
        var grp = try Kn.Group.init(gpa, ps[0..count]);
        hits.clearRetainingCapacity();
        try grp.scan(gpa, g, max_cost, &hits);
        for (hits.items) |h| {
            if (h.lane < count) by_guide[start + h.lane] += 1;
        }
        total_hits += hits.items.len;
        grp.deinit(gpa);
    }
    const sec = @as(f64, @floatFromInt(now(io) - t0)) * 1e-9;
    var reported: usize = 0;
    for (by_guide) |x| reported += x;
    const cells = @as(f64, @floatFromInt(g.len)) * guide_len * @as(f64, @floatFromInt(nguides));
    std.debug.print("genome {d} bp, {d} guides x {d} nt, K={d}, lanes={d}{s}: {d} hits in {d:.3} s, {d:.2} GCUPS\n", .{
        g.len, nguides, guide_len, max_cost, Kn.lanes, if (pack) " packed" else "", reported, sec, cells / sec * 1e-9,
    });

    if (!verify) return;
    // Scalar check: the last DP row at each column is the best cost ending there.
    const t1 = now(io);
    var expected: usize = 0;
    var mismatched: usize = 0;
    for (guides, 0..) |gd, gi| {
        var col: [guide_len + 1]i64 = undefined;
        for (&col, 0..) |*x, i| x.* = @as(i64, @intCast(i)) * scheme.gap;
        var n: usize = 0;
        for (g) |c| {
            var diag = col[0];
            col[0] = 0;
            for (gd, 1..) |pc, i| {
                const cell = @min(diag + scheme.cost(pc, c), @min(col[i] + scheme.gap, col[i - 1] + scheme.gap));
                diag = col[i];
                col[i] = cell;
            }
            if (col[guide_len] <= max_cost) n += 1;
        }
        expected += n;
        if (n != by_guide[gi]) mismatched += 1;
    }
    const sec_scalar = @as(f64, @floatFromInt(now(io) - t1)) * 1e-9;
    std.debug.print("scalar check: {d} hits expected, {d} guides disagree; scalar DP took {d:.3} s ({d:.2} GCUPS)\n", .{
        expected, mismatched, sec_scalar, cells / sec_scalar * 1e-9,
    });
    if (mismatched != 0 or expected != reported) return error.Mismatch;
}
