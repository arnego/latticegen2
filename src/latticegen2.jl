module latticegen2

# Module layout follows docs/algorithm.md §12. Include order matters: later
# files depend on error types / helpers defined in earlier ones.
include("runlog.jl")     # error taxonomy, exit codes, logging
include("cli.jl")        # argument parsing/validation
include("lattice.jl")    # lattice math: directions, basis, strut enumeration
include("geomkernel.jl") # gmsh/OCCT primitives: I/O, prototypes, booleans
include("classify.jl")   # boundary classification: mesh, spatial hash, distance, ray-cast
include("tiling.jl")     # tile partitioning and periodicity-shortcut detection
include("pipeline.jl")   # top-level orchestration, calibration, Distributed workers
include("stepmeta.jl")   # STEP header rewrite, round-trip verification
include("app_entry.jl")  # PackageCompiler.jl standalone-app entry point (specification.md §2)

export
    # app_entry.jl
    julia_main,
    # runlog.jl
    LatticeGenError, ParamError, InputGeometryError, ProcessingError,
    ResourceError, OutputError, CancelledError, exit_code, is_interrupt,
    RunLog, open_runlog, close_runlog, log_line, log_always, log_stage,
    log_run_header, log_summary, log_failure, log_cancelled, update_rss!, observe_rss!,
    format_param, format_bytes, format_duration,
    # cli.jl
    CliArgs, parse_args, resolve_output_paths, preflight_checks, USAGE,
    # lattice.jl
    Vec3, dot3, cross3, norm3, normalize3,
    THETA, strut_direction, cube_edge, LatticeParams, node, profile_vertices,
    StrutRef, endpoints, midpoint, index_range, candidate_struts,
    aabbs_overlap, overlap_components,
    # geomkernel.jl
    with_gmsh, capture_output, import_shapes, bounding_box, write_model,
    sync_model, bounding_box_nosync, bounding_boxes, remove_model_entities, model_solids,
    new_model, build_prototype, build_prototypes, instantiate_strut,
    fuse_all, common_with, cut_with, volume_of, remove_entities, copy_translate,
    # classify.jl
    StrutClass, INTERIOR, BOUNDARY, OUTSIDE,
    TriMesh, mesh_chordal_target, tessellate_surface, tri_aabb,
    SpatialHash, build_spatial_hash, query_cells,
    closest_point_segment, dist_point_segment, closest_point_triangle,
    dist_point_triangle, dist_segment_segment, dist_segment_triangle,
    segment_mesh_min_distance, segment_triangle_intersect,
    ray_triangle, raycast_parity, mesh_bounds, point_inside, classify_strut,
    # tiling.jl
    TileKey, tile_key, TileDescriptor, partition_tiles, is_full_interior,
    tile_translation, interface_nodes,
    # pipeline.jl
    make_run_tempdir, set_below_normal_priority!, shutdown_workers!, tile_brep_path, part_name,
    determine_worker_count, calibrate, derive_tile_size, process_tile, trim_disjoint,
    drop_tile_local_islands, DEFAULT_TILE_FUSE_BUDGET_SECONDS, MEMORY_WATCHDOG_MAX_PAUSE_SECONDS,
    dispatch_jobs, dispatch_tiles, balanced_fuse!, filter_floating!, assert_no_stray_solids,
    merge_group, do_merge_round, assemble!, run_pipeline,
    # stepmeta.jl
    generation_params_text, rewrite_step_header!, round_trip_check

end # module latticegen2
