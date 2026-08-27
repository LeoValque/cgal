// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SOUP_AND_TREE_HELPER_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SOUP_AND_TREE_HELPER_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/property_map.h>
#include <fstream>
#include <sstream>
#include <set>
#include <type_traits>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<class TriangleMesh, class GT>
struct AABB_tree_build_helper
{
  using Primitive = AABB_face_graph_triangle_primitive<TriangleMesh>;
  using Traits = AABB_traits_3<GT, Primitive>;
  using Tree = AABB_tree<Traits>;

  using Graph_traits = boost::graph_traits<TriangleMesh>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;

  template <class RPM>
  struct Split_primitives
  {
    Split_primitives(const RPM &rpm): rpm(rpm){}

    template<typename PrimitiveIterator>
    void operator()(PrimitiveIterator first,
                    PrimitiveIterator beyond,
                    const CGAL::Bbox_3& bbox) const
      {
        auto longest_axis=[](const CGAL::Bbox_3& bbox){
          const double dx = bbox.x_span();
          const double dy = bbox.y_span();
          const double dz = bbox.z_span();
          return (dx>=dy) ? ((dx>=dz) ? 0 : 2) : ((dy>=dz) ? 1 : 2);
        };

        PrimitiveIterator middle = first + (beyond - first)/2;
        const int crd=longest_axis(bbox);
        std::nth_element(first, middle, beyond, [this, crd](const Primitive& p1, const Primitive& p2){ return get(rpm, p1.id())[crd] < get(rpm, p2.id())[crd];});
      }
    const RPM &rpm;
  };

  // For exact side_of_triangle_mesh
  template <class BPM>
  struct Compute_bbox {
    Compute_bbox(const BPM& bpm): bpm(bpm){}

    template<typename ConstPrimitiveIterator>
    CGAL::Bbox_3 operator()(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator beyond) const
    {
      CGAL::Bbox_3 bbox = get(bpm, first->id());
      for(++first; first != beyond; ++first)
        bbox += get(bpm, first->id());
      return bbox;
    }
    BPM bpm;
  };

  template<class Concurrency_tag=Sequential_tag, class FaceRange, class VertexPointMap>
  void build(const FaceRange &face_range, Tree& tree, const TriangleMesh& tm, VertexPointMap& vpm){
    using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>;
    using Face_ref_point_tag = typename CGAL::dynamic_face_property_t<Epick::Point_3>;
    using Bbox_map = typename boost::property_map<TriangleMesh, Face_bbox_tag>::const_type;
    using Ref_point_map = typename boost::property_map<TriangleMesh, Face_ref_point_tag>::const_type;

    using VPM_kernel = typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel;
    CGAL::Cartesian_converter<VPM_kernel, Epick> to_input;

    Bbox_map bb_map = get(Face_bbox_tag(), tm);
    Ref_point_map rp_map = get(Face_ref_point_tag(), tm);

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
    {
      tbb::parallel_for(std::size_t(0), faces(tm).size(), [&](std::size_t i){
        face_descriptor f(i);
        put(bb_map, f, face_bbox(f, tm));
        put(rp_map, f, to_input(get(vpm, target(halfedge(f, tm), tm))) );
      });
    }
    else
#endif
    {
      for(face_descriptor f : faces(tm)){
        put(bb_map, f, face_bbox(f, tm));
        put(rp_map, f, to_input(get(vpm, target(halfedge(f, tm), tm))) );
      }
    }
    tree.insert(face_range.begin(), face_range.end(), tm, vpm);

    Compute_bbox<Bbox_map> compute_bbox(bb_map);
    Split_primitives<Ref_point_map> split_primitives(rp_map);
    tree.template custom_build<Concurrency_tag>(compute_bbox, split_primitives);
  }

  template<class ConcurrencyTag=Sequential_tag, class VertexPointMap>
  void build(Tree& tree, const TriangleMesh& tm, VertexPointMap& vpm){
    build<ConcurrencyTag>(faces(tm), tree, tm, vpm);
  }
};

} } } // CGAL::Polygon_mesh_processing::internal

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SOUP_AND_TREE_HELPER_H
