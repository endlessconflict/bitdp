//! Ablation: word operations per column for one scheme under each builder
//! configuration. One scheme per build (compile memory grows with every
//! derivation): zig build ablation -Dscheme=N, N indexes `list` below.

const std = @import("std");
const bitdp = @import("bitdp");
const options = @import("options");

const list = [_]struct { []const u8, bitdp.Scheme }{
    .{ "edit (0,1,1)", bitdp.schemes.edit },
    .{ "indel (0,2,1)", .{ .mismatch = 2 } },
    .{ "(0,1,2)", .{ .mismatch = 1, .gap = 2 } },
    .{ "(0,3,2)", .{ .mismatch = 3, .gap = 2 } },
    .{ "bitpal (2,-3,-5)", bitdp.schemes.bitpal[1] },
    .{ "bitpal (3,-4,-6)", bitdp.schemes.bitpal[2] },
    .{ "bitpal (4,-5,-9)", bitdp.schemes.bitpal[3] },
    .{ "bitpal (4,-7,-11)", bitdp.schemes.bitpal[4] },
    .{ "tstv", bitdp.schemes.ts_tv },
    .{ "blosum62, gap 4", bitdp.schemes.blosum62Linear(4) },
};

const configs = [_]struct { []const u8, bitdp.Options }{
    .{ "all", .{} },
    .{ "no merge networks", .{ .merge = false } },
    .{ "no truth tables", .{ .truth_table = false } },
    .{ "direct only", .{ .truth_table = false, .merge = false } },
    // Bergeron and Hamel's algorithm on bit planes, compiled by the same back end.
    .{ "B&H planes", bitdp.Options.bh_only },
    .{ "best of both", bitdp.Options.best_of_both },
};

/// Scaling sweep, -Dscheme=100+i: costs (0, 1, g) keep the chained block
/// small while the difference range grows; (-1, g, g+1) grow both, like
/// BitPAl's weight sets.
const sweep = blk: {
    var s: [24]struct { []const u8, bitdp.Scheme } = undefined;
    for (1..13) |g| {
        s[2 * (g - 1)] = .{ std.fmt.comptimePrint("sweep (0,1,{d})", .{g}), .{ .mismatch = 1, .gap = g } };
        s[2 * (g - 1) + 1] = .{ std.fmt.comptimePrint("sweep (-1,{d},{d})", .{ g, g + 1 }), .{ .match = -1, .mismatch = g, .gap = g + 1 } };
    }
    break :blk s;
};

/// Nodes per op kind among those the outputs use (the nodes Table 1 counts),
/// so the per-word cost of the multi-word path can be recomputed from it.
fn opMix(plan: bitdp.Plan) [@typeInfo(bitdp.Op).@"enum".fields.len]u32 {
    @setEvalBranchQuota(1 << 20);
    var live = [_]bool{false} ** plan.nodes.len;
    for (plan.out_v[0..plan.nstate]) |o| live[o] = true;
    for (plan.out_u[0..plan.nstate]) |o| live[o] = true;
    var mix = [_]u32{0} ** @typeInfo(bitdp.Op).@"enum".fields.len;
    var i: usize = plan.len;
    while (i > 0) {
        i -= 1;
        if (!live[i]) continue;
        const n = plan.nodes[i];
        mix[@intFromEnum(n.op)] += 1;
        switch (n.op) {
            .input, .zero, .ones => {},
            .not, .shl0, .shl1 => live[n.a] = true,
            else => {
                live[n.a] = true;
                live[n.b] = true;
            },
        }
    }
    return mix;
}

pub fn main() void {
    if (options.scheme >= 100) {
        const entry = sweep[if (options.scheme >= 100) options.scheme - 100 else 0];
        const plan = comptime bitdp.derive(entry[1]);
        const bh = comptime bitdp.deriveWith(entry[1], bitdp.Options.bh_only);
        std.debug.print("{s}: k={d} chains={d} ops={d} bh={d}\n", .{ entry[0], plan.k, plan.chains, plan.cost, bh.cost });
        return;
    }
    const entry = list[if (options.scheme >= 100) 0 else options.scheme];
    const mix = comptime opMix(bitdp.derive(entry[1]));
    std.debug.print("{s} live ops:", .{entry[0]});
    inline for (@typeInfo(bitdp.Op).@"enum".fields, 0..) |f, i| {
        if (mix[i] > 0) std.debug.print(" {s}={d}", .{ f.name, mix[i] });
    }
    std.debug.print("\n{s}:", .{entry[0]});
    inline for (configs) |c| {
        const opt = c[1];
        const usable = opt.direct or opt.merge or opt.truth_table or opt.bh;
        if (usable) {
            const plan = comptime bitdp.deriveWith(entry[1], opt);
            std.debug.print("  {s}={d}", .{ c[0], plan.cost });
        }
    }
    // Their listing, every printed operation executed (the count review/factcheck/bh.py reports).
    std.debug.print("  B&H listing={d}\n", .{comptime bitdp.bhListing(entry[1]).cost});
}
