//! C interface: precompiled kernels for a fixed set of schemes (include/bitdp.h).
//! Custom schemes stay Zig-only, because the derivation runs at compile time.
//!
//! Every function returns 0 on success, BITDP_EINVAL (1) for a bad argument and
//! BITDP_ENOMEM (2) when an allocation fails; results go through pointers. Memory the
//! library hands out is freed with the library's own *_free functions.

const std = @import("std");
const bitdp = @import("root.zig");

const gpa = std.heap.smp_allocator;
const ok: c_int = 0;
const einval: c_int = 1;
const enomem: c_int = 2;

/// A search-mode hit: pattern `pattern` (index into the caller's list) ends at text
/// position `end` (exclusive, 1-based) with alignment cost `cost`.
pub const Hit = extern struct { pattern: usize, end: usize, cost: i64 };

fn bytes(p: ?[*]const u8, n: usize) ?[]const u8 {
    if (n == 0) return &.{};
    return (p orelse return null)[0..n];
}

/// Slices for a C list of (pointer, length) pairs; null when any entry is invalid.
fn slices(ps: ?[*]const ?[*]const u8, lens: ?[*]const usize, count: usize, min_len: usize) error{ OutOfMemory, Invalid }![][]const u8 {
    const p = ps orelse return error.Invalid;
    const l = lens orelse return error.Invalid;
    const out = try gpa.alloc([]const u8, count);
    errdefer gpa.free(out);
    for (out, 0..) |*s, i| {
        if (l[i] < min_len) return error.Invalid;
        s.* = bytes(p[i], l[i]) orelse return error.Invalid;
    }
    return out;
}

fn Api(comptime scheme: bitdp.Scheme) type {
    const K = bitdp.Kernel(scheme);
    return struct {
        fn distance(p: ?[*]const u8, m: usize, t: ?[*]const u8, n: usize, out: ?*i64) callconv(.c) c_int {
            if (m < 1 or m > 62) return einval;
            const pat = bytes(p, m) orelse return einval;
            const text = bytes(t, n) orelse return einval;
            (out orelse return einval).* = K.distance(pat, text);
            return ok;
        }

        fn alignerNew(p: ?[*]const u8, m: usize, out: ?*?*K.Aligner) callconv(.c) c_int {
            const dst = out orelse return einval;
            if (m < 1) return einval;
            const pat = bytes(p, m) orelse return einval;
            const a = gpa.create(K.Aligner) catch return enomem;
            a.* = K.Aligner.init(gpa, pat) catch {
                gpa.destroy(a);
                return enomem;
            };
            dst.* = a;
            return ok;
        }

        fn alignerDistance(a: ?*K.Aligner, t: ?[*]const u8, n: usize, out: ?*i64) callconv(.c) c_int {
            const al = a orelse return einval;
            const text = bytes(t, n) orelse return einval;
            (out orelse return einval).* = al.distance(text);
            return ok;
        }

        fn alignerFree(a: ?*K.Aligner) callconv(.c) void {
            const al = a orelse return;
            al.deinit(gpa);
            gpa.destroy(al);
        }

        fn distances(ps: ?[*]const ?[*]const u8, ms: ?[*]const usize, ts: ?[*]const ?[*]const u8, ns: ?[*]const usize, count: usize, out: ?[*]i64) callconv(.c) c_int {
            if (count == 0) return ok;
            const o = out orelse return einval;
            const pats = slices(ps, ms, count, 1) catch |e| return if (e == error.OutOfMemory) enomem else einval;
            defer gpa.free(pats);
            const texts = slices(ts, ns, count, 0) catch |e| return if (e == error.OutOfMemory) enomem else einval;
            defer gpa.free(texts);
            K.distances(gpa, pats, texts, o[0..count]) catch return enomem;
            return ok;
        }

        fn scan(ps: ?[*]const ?[*]const u8, ms: ?[*]const usize, count: usize, t: ?[*]const u8, n: usize, max_cost: i64, hits_out: ?*?[*]Hit, nhits: ?*usize) callconv(.c) c_int {
            const ho = hits_out orelse return einval;
            const no = nhits orelse return einval;
            ho.* = null;
            no.* = 0;
            if (count == 0) return ok;
            const text = bytes(t, n) orelse return einval;
            const pats = slices(ps, ms, count, 1) catch |e| return if (e == error.OutOfMemory) enomem else einval;
            defer gpa.free(pats);
            var all: std.ArrayList(Hit) = .empty;
            defer all.deinit(gpa);
            var hits: std.ArrayList(K.Hit) = .empty;
            defer hits.deinit(gpa);
            var start: usize = 0;
            while (start < count) : (start += K.lanes) {
                const k = @min(K.lanes, count - start);
                var g = K.Group.init(gpa, pats[start..][0..k]) catch return enomem;
                defer g.deinit(gpa);
                hits.clearRetainingCapacity();
                g.scan(gpa, text, max_cost, &hits) catch return enomem;
                for (hits.items) |h| {
                    if (h.lane < k) all.append(gpa, .{ .pattern = start + h.lane, .end = h.end, .cost = h.cost }) catch return enomem;
                }
            }
            if (all.items.len == 0) return ok;
            const owned = all.toOwnedSlice(gpa) catch return enomem;
            ho.* = owned.ptr;
            no.* = owned.len;
            return ok;
        }
    };
}

fn hitsFree(hits: ?[*]Hit, n: usize) callconv(.c) void {
    const h = hits orelse return;
    gpa.free(h[0..n]);
}

fn searchMode(comptime s: bitdp.Scheme) bitdp.Scheme {
    var t = s;
    t.mode = .search;
    return t;
}

const exported = .{
    .{ "edit_global", bitdp.schemes.edit },
    .{ "edit_search", searchMode(bitdp.schemes.edit) },
    .{ "tstv_global", bitdp.schemes.ts_tv },
    .{ "tstv_search", searchMode(bitdp.schemes.ts_tv) },
    .{ "bitpal_global", bitdp.schemes.bitpal[1] },
};

comptime {
    for (exported) |e| {
        const A = Api(e[1]);
        const pre = "bitdp_" ++ e[0];
        @export(&A.distance, .{ .name = pre ++ "_distance" });
        @export(&A.alignerNew, .{ .name = pre ++ "_aligner_new" });
        @export(&A.alignerDistance, .{ .name = pre ++ "_aligner_distance" });
        @export(&A.alignerFree, .{ .name = pre ++ "_aligner_free" });
        @export(&A.distances, .{ .name = pre ++ "_distances" });
        if (e[1].mode == .search) @export(&A.scan, .{ .name = pre ++ "_scan" });
    }
    @export(&hitsFree, .{ .name = "bitdp_hits_free" });
}
