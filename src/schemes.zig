//! Ready-made scoring schemes, as costs to minimize.

const Scheme = @import("derive.zig").Scheme;

/// Unit edit distance.
pub const edit: Scheme = .{};

/// The weight sets BitPAl benchmarks, score (M, I, G) written as costs
/// (-M, -I, -G), which gives the same optimal alignments.
pub const bitpal = [_]Scheme{
    .{ .match = 0, .mismatch = 1, .gap = 1 },
    .{ .match = -2, .mismatch = 3, .gap = 5 },
    .{ .match = -3, .mismatch = 4, .gap = 6 },
    .{ .match = -4, .mismatch = 5, .gap = 9 },
    .{ .match = -4, .mismatch = 7, .gap = 11 },
};

/// DNA with transitions (A<->G, C<->T) cheaper than transversions.
pub fn tsTvCost(a: u8, b: u8) i32 {
    if (a == b) return 0;
    const pa = a == 'A' or a == 'G';
    const pb = b == 'A' or b == 'G';
    return if (pa == pb) 1 else 2;
}

pub const ts_tv: Scheme = .{ .sub = &tsTvCost, .gap = 2 };

/// BLOSUM62 scores (NCBI, ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62), 20 standard residues.
pub const blosum62_order = "ARNDCQEGHILKMFPSTWYV";
pub const blosum62 = [20][20]i8{
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
pub fn blosum62Cost(a: u8, b: u8) i32 {
    return -@as(i32, blosum62[blosum_index[a]][blosum_index[b]]);
}

/// BLOSUM62 with a linear gap cost.
pub fn blosum62Linear(gap: i32) Scheme {
    return .{ .sub = &blosum62Cost, .alphabet = blosum62_order, .gap = gap };
}
