

def assert_view_contains_pharmacophore(view, n_points):
    # The protein is there
    assert "nglview.adaptor.MDTrajTrajectory" in view._ngl_component_names
    # And the ligand as well
    assert "nglview.adaptor.RdkitStructure" in view._ngl_component_names
    # There is at least one sphere
    assert "nglview.shape.Shape" in view._ngl_component_names

    n_shapes = len([comp for comp in view._ngl_component_names
                    if comp == "nglview.shape.Shape"])
    assert n_shapes >= n_points  # Direction arrows can make n_shapes greater
