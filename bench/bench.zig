//! bitdp side of the benchmark: global alignment cost of every pair in a
//! "pattern<TAB>text" file. Output matches bench/baselines.c.
//!
//! usage: bitdp-bench {edit|bitpal|tstv|blosum} FILE [batch]
//!
//! With "batch", pairs go through Kernel.distances, one alignment per SIMD lane.

const std = @import("std");
const bitdp = @import("bitdp");

fn runBatch(comptime s: bitdp.Scheme, name: []const u8, io: std.Io, gpa: std.mem.Allocator, data: []const u8) !void {
    const K = bitdp.Kernel(s);
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
    const t0 = std.Io.Timestamp.now(io, .awake).nanoseconds;
    try K.distances(gpa, ps.items, ts.items, out);
    const sec = @as(f64, @floatFromInt(std.Io.Timestamp.now(io, .awake).nanoseconds - t0)) * 1e-9;
    var total: i64 = 0;
    for (out) |x| total += x;
    std.debug.print("bitdp-batch{d} {s} pairs={d} cost_sum={d} seconds={d:.4} gcups={d:.3} ops_per_word={d}\n", .{
        K.lanes, name, out.len, total, sec, @as(f64, @floatFromInt(cells)) / sec * 1e-9, K.ops,
    });
}

fn run(comptime s: bitdp.Scheme, name: []const u8, io: std.Io, data: []const u8) !void {
    const K = bitdp.Kernel(s);
    var arena_buf: [1 << 20]u8 = undefined;
    var total: i64 = 0;
    var cells: u64 = 0;
    var pairs: u64 = 0;
    var busy: i96 = 0;
    var lines = std.mem.tokenizeScalar(u8, data, '\n');
    while (lines.next()) |line| {
        const tab = std.mem.indexOfScalar(u8, line, '\t') orelse continue;
        const p = line[0..tab];
        const t = line[tab + 1 ..];
        var fba = std.heap.FixedBufferAllocator.init(&arena_buf);
        const t0 = std.Io.Timestamp.now(io, .awake).nanoseconds;
        var a = try K.Aligner.init(fba.allocator(), p);
        total += a.distance(t);
        busy += std.Io.Timestamp.now(io, .awake).nanoseconds - t0;
        cells += p.len * t.len;
        pairs += 1;
    }
    const sec = @as(f64, @floatFromInt(busy)) * 1e-9;
    std.debug.print("bitdp {s} pairs={d} cost_sum={d} seconds={d:.4} gcups={d:.3} ops_per_word={d}\n", .{
        name, pairs, total, sec, @as(f64, @floatFromInt(cells)) / sec * 1e-9, K.ops,
    });
}

const table = .{
    .{ "edit", bitdp.schemes.edit },
    .{ "bitpal", bitdp.schemes.bitpal[1] },
    .{ "tstv", bitdp.schemes.ts_tv },
    .{ "blosum", bitdp.schemes.blosum62Linear(4) },
};

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    if (args.len < 3) {
        std.debug.print("usage: bitdp-bench SCHEME FILE\n", .{});
        return error.Usage;
    }
    const data = try std.Io.Dir.cwd().readFileAlloc(io, args[2], arena, .unlimited);
    const scheme = args[1];
    inline for (table) |e| {
        if (std.mem.eql(u8, scheme, e[0])) {
            if (args.len > 3) return runBatch(e[1], e[0], io, arena, data);
            return run(e[1], e[0], io, data);
        }
    }
    return error.UnknownScheme;
}
