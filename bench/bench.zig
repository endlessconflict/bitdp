//! bitdp side of the benchmark: global alignment cost of every pair in a
//! "pattern<TAB>text" file. Output matches bench/baselines.c, including the
//! per-pair cost hash of bench/pairhash.h and the BENCH_DUMP file.
//!
//! usage: bitdp-bench {edit|bitpal|tstv|blosum} FILE [batch]
//!
//! With "batch", pairs go through Kernel.distances, one alignment per SIMD lane.
//!
//! usage: bitdp-bench SCHEME ava QUERIES SUBJECTS
//!
//! All queries against all subjects, one sequence per line (BGSA's setting):
//! subjects are prepared once in lane groups, queries are the shared texts.
//!
//! A scheme name with the prefix "bh-" (bh-edit, ...) runs the same kernel on
//! Bergeron and Hamel's algorithm instead of the derived program.

const std = @import("std");
const bitdp = @import("bitdp");

/// Prints the run's line: cost sum and FNV-1a (64 bit) over each cost in decimal plus a newline, which
/// is also the content of the file BENCH_DUMP names, if set. Runs after the timed region.
fn report(io: std.Io, gpa: std.mem.Allocator, env: *std.process.Environ.Map, label: []const u8, lanes: usize, name: []const u8, costs: []const i64, cells: u64, sec: f64, rounds: u64, ops: usize) !void {
    var text: std.ArrayList(u8) = .empty;
    var h: u64 = 0xcbf29ce484222325;
    var total: i64 = 0;
    for (costs) |c| {
        const start = text.items.len;
        try text.print(gpa, "{d}\n", .{c});
        for (text.items[start..]) |b| h = (h ^ b) *% 0x100000001b3;
        total += c;
    }
    if (env.get("BENCH_DUMP")) |path| try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = path, .data = text.items });
    std.debug.print("{s}{d} {s} pairs={d} cost_sum={d} hash={x:0>16} seconds={d:.4} gcups={d:.3} ops_per_word={d} rounds={d}\n", .{
        label, lanes, name, costs.len, total, h, sec, @as(f64, @floatFromInt(cells * rounds)) / sec * 1e-9, ops, rounds,
    });
}

/// The timed region of every mode is the whole loop over the input, repeated until BENCH_MIN_SECONDS
/// (default 0.25) have passed; the rate is over all rounds (the same rule as bench/pairhash.h).
fn minSeconds(env: *std.process.Environ.Map) f64 {
    return std.fmt.parseFloat(f64, env.get("BENCH_MIN_SECONDS") orelse "0.25") catch 0.25;
}

fn since(io: std.Io, t0: i96) f64 {
    return @as(f64, @floatFromInt(std.Io.Timestamp.now(io, .awake).nanoseconds - t0)) * 1e-9;
}

fn runBatch(comptime s: bitdp.Scheme, comptime opt: bitdp.Options, name: []const u8, io: std.Io, gpa: std.mem.Allocator, env: *std.process.Environ.Map, data: []const u8) !void {
    const K = bitdp.KernelWith(s, opt);
    var ps: std.ArrayList([]const u8) = .empty;
    var ts: std.ArrayList([]const u8) = .empty;
    var cells: u64 = 0;
    var lines = std.mem.tokenizeScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const tab = std.mem.indexOfScalar(u8, line, '\t') orelse continue;
        try ps.append(gpa, line[0..tab]);
        try ts.append(gpa, line[tab + 1 ..]);
        cells += tab * (line.len - tab - 1);
    }
    const out = try gpa.alloc(i64, ps.items.len);
    const min = minSeconds(env);
    var rounds: u64 = 0;
    var sec: f64 = 0;
    const t0 = std.Io.Timestamp.now(io, .awake).nanoseconds;
    while (rounds == 0 or sec < min) : (rounds += 1) {
        try K.distances(gpa, ps.items, ts.items, out);
        sec = since(io, t0);
    }
    try report(io, gpa, env, "bitdp-batch", K.lanes, name, out, cells, sec, rounds, K.ops);
}

fn runAllVsAll(comptime s: bitdp.Scheme, comptime opt: bitdp.Options, name: []const u8, io: std.Io, gpa: std.mem.Allocator, env: *std.process.Environ.Map, qdata: []const u8, sdata: []const u8) !void {
    const K = bitdp.KernelWith(s, opt);
    var queries: std.ArrayList([]const u8) = .empty;
    var subjects: std.ArrayList([]const u8) = .empty;
    var it = std.mem.tokenizeScalar(u8, qdata, '\n');
    while (it.next()) |l| try queries.append(gpa, l);
    it = std.mem.tokenizeScalar(u8, sdata, '\n');
    while (it.next()) |l| try subjects.append(gpa, l);
    const ns = subjects.items.len;
    const costs = try gpa.alloc(i64, queries.items.len * ns); // query-major

    const min = minSeconds(env);
    var rounds: u64 = 0;
    var sec: f64 = 0;
    const t0 = std.Io.Timestamp.now(io, .awake).nanoseconds;
    while (rounds == 0 or sec < min) : (rounds += 1) { // a round prepares the subjects again
        const ngroups = (ns + K.lanes - 1) / K.lanes;
        const groups = try gpa.alloc(K.Group, ngroups);
        for (groups, 0..) |*g, i| {
            const lo = i * K.lanes;
            g.* = try K.Group.init(gpa, subjects.items[lo..@min(lo + K.lanes, ns)]);
        }
        for (queries.items, 0..) |q, qi| for (groups, 0..) |*g, i| {
            const r = g.distances(q);
            const count = @min(K.lanes, ns - i * K.lanes);
            @memcpy(costs[qi * ns + i * K.lanes ..][0..count], r[0..count]);
        };
        sec = since(io, t0);
    }
    var slen: u64 = 0;
    for (subjects.items) |x| slen += x.len;
    var cells: u64 = 0;
    for (queries.items) |q| cells += q.len * slen;
    try report(io, gpa, env, "bitdp-group", K.lanes, name, costs, cells, sec, rounds, K.ops);
}

fn run(comptime s: bitdp.Scheme, comptime opt: bitdp.Options, name: []const u8, io: std.Io, gpa: std.mem.Allocator, env: *std.process.Environ.Map, data: []const u8) !void {
    const K = bitdp.KernelWith(s, opt);
    var arena_buf: [1 << 20]u8 = undefined;
    var ps: std.ArrayList([]const u8) = .empty;
    var ts: std.ArrayList([]const u8) = .empty;
    var cells: u64 = 0;
    var lines = std.mem.tokenizeScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const tab = std.mem.indexOfScalar(u8, line, '\t') orelse continue;
        try ps.append(gpa, line[0..tab]);
        try ts.append(gpa, line[tab + 1 ..]);
        cells += tab * (line.len - tab - 1);
    }
    const costs = try gpa.alloc(i64, ps.items.len);
    const min = minSeconds(env);
    var rounds: u64 = 0;
    var sec: f64 = 0;
    const t0 = std.Io.Timestamp.now(io, .awake).nanoseconds;
    while (rounds == 0 or sec < min) : (rounds += 1) {
        for (ps.items, ts.items, costs) |p, t, *c| {
            var fba = std.heap.FixedBufferAllocator.init(&arena_buf);
            var a = try K.Aligner.init(fba.allocator(), p);
            c.* = a.distance(t);
        }
        sec = since(io, t0);
    }
    try report(io, gpa, env, "bitdp", 1, name, costs, cells, sec, rounds, K.ops);
}

const table = .{
    .{ "edit", bitdp.schemes.edit },
    .{ "bitpal", bitdp.schemes.bitpal[1] },
    .{ "tstv", bitdp.schemes.ts_tv },
    .{ "blosum", bitdp.schemes.blosum62Linear(4) },
};

/// The derived program, and Bergeron and Hamel's algorithm for the "bh-" names.
const builders = [_]bitdp.Options{ .{}, bitdp.Options.bh_only };

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const arena = init.arena.allocator();
    const env = init.environ_map;
    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 3) {
        std.debug.print("usage: bitdp-bench SCHEME FILE\n", .{});
        return error.Usage;
    }
    const bh = std.mem.startsWith(u8, args[1], "bh-");
    const scheme = args[1][if (bh) 3 else 0..];
    if (std.mem.eql(u8, args[2], "ava")) {
        if (args.len != 5) return error.Usage;
        const q = try std.Io.Dir.cwd().readFileAlloc(io, args[3], arena, .unlimited);
        const sb = try std.Io.Dir.cwd().readFileAlloc(io, args[4], arena, .unlimited);
        inline for (table) |e| inline for (builders) |b| {
            if (bh == b.bh and std.mem.eql(u8, scheme, e[0])) return runAllVsAll(e[1], b, args[1], io, arena, env, q, sb);
        };
        return error.UnknownScheme;
    }
    const data = try std.Io.Dir.cwd().readFileAlloc(io, args[2], arena, .unlimited);
    inline for (table) |e| inline for (builders) |b| {
        if (bh == b.bh and std.mem.eql(u8, scheme, e[0])) {
            if (args.len > 3) return runBatch(e[1], b, args[1], io, arena, env, data);
            return run(e[1], b, args[1], io, arena, env, data);
        }
    };
    return error.UnknownScheme;
}
