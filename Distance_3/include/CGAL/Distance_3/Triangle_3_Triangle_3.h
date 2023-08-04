// Copyright (c) 1998-2021
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_DISTANCE_3_TRIANGLE_3_TRIANGLE_3_H
#define CGAL_DISTANCE_3_TRIANGLE_3_TRIANGLE_3_H

#include <CGAL/Distance_3/Point_3_Point_3.h>
#include <CGAL/Distance_3/Segment_3_Segment_3.h>
#include <CGAL/Intersections_3/Triangle_3_Triangle_3.h>

#include <CGAL/Triangle_3.h>

namespace CGAL {
namespace Distance_3 {
namespace internal {

template <typename K>
std::pair<Segment_3_Segment_3_Result<K>, bool>
test_edge_pair(const typename K::Point_3& p1,
               const typename K::Point_3& q1,
               const typename K::Point_3& r1,
               const typename K::Point_3& p2,
               const typename K::Point_3& q2,
               const typename K::Point_3& r2,
               const K& k,
               typename K::FT& global_min_sqd,
               bool& are_triangles_known_to_be_disjoint)
{
  typedef typename K::FT                                                  FT;
  typedef typename K::Point_3                                             Point_3;
  typedef typename K::Vector_3                                            Vector_3;

  typename K::Compute_scalar_product_3 scalar_product = k.compute_scalar_product_3_object();
  typename K::Construct_segment_3 segment = k.construct_segment_3_object();
  typename K::Construct_scaled_vector_3 scale_vector = k.construct_scaled_vector_3_object();
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Construct_translated_point_3 translate = k.construct_translated_point_3_object();

  Distance_3::internal::Segment_3_Segment_3_Result<K> res =
      internal::squared_distance(segment(p1, q1), segment(p2, q2), k);

  if(res.squared_distance <= global_min_sqd)
    global_min_sqd = res.squared_distance;
  else
    return std::make_pair(res, false);

  const Vector_3 v1 = vector(p1, q1), v2 = vector(p2, q2);
  const Point_3 m1 = translate(p1, scale_vector(v1, res.x));
  const Point_3 m2 = translate(p2, scale_vector(v2, res.y));
  const Vector_3 vr1 = vector(m1, r1), vr2 = vector(m2, r2);
  const Vector_3 n = vector(m1, m2);

  const FT sp_r1 = scalar_product(vr1, n);
  const FT sp_r2 = scalar_product(vr2, n);
  const bool is_r1_closer = (sp_r1 > 0); // Plane_3{m1, n}.has_on_positive_side(r1);
  const bool is_r2_closer = (sp_r2 < 0); // Plane_3{m2, -n}.has_on_positive_side(r2);
  const bool is_best_pair = !is_r1_closer && !is_r2_closer;

  // Even if it is not the best pair, one may be able to deduce if the triangles do not intersect
  // by checking if there is a void space between the planes orthogonal to the vector realizing
  // the min distance between the edges and passing through the third points.
  if(!is_best_pair)
  {
    FT separating_distance = res.squared_distance;
    if(is_r1_closer)
      separating_distance -= sp_r1;
    if(is_r2_closer)
      separating_distance += sp_r2;

    if(separating_distance > 0)
      are_triangles_known_to_be_disjoint = true;
  }

  return std::make_pair(res, is_best_pair);
}

template <typename K>
std::pair<typename K::FT, bool>
test_vertex_triangle(const typename K::Triangle_3& tr1,
                     const typename K::Triangle_3& tr2,
                     const K& k,
                     bool& are_triangles_known_to_be_disjoint)
{
  typedef typename K::FT                                                  FT;
  typedef typename K::Point_3                                             Point_3;
  typedef typename K::Vector_3                                            Vector_3;

  typename K::Compute_scalar_product_3 scalar_product = k.compute_scalar_product_3_object();
  typename K::Construct_cross_product_vector_3 cross_product = k.construct_cross_product_vector_3_object();
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();

  const Point_3& p1 = vertex(tr1, 0);
  const Point_3& q1 = vertex(tr1, 1);
  const Point_3& r1 = vertex(tr1, 2);
  const Point_3& p2 = vertex(tr2, 0);
  const Point_3& q2 = vertex(tr2, 1);
  const Point_3& r2 = vertex(tr2, 2);

  const Vector_3 p2q2 = vector(p2, q2);
  const Vector_3 p2r2 = vector(p2, r2);
  const Vector_3 n2 = cross_product(p2q2, p2r2);

  if(scalar_product(n2, n2) == FT(0))
    return std::make_pair(0, false);

  std::array<FT, 3> sps = { scalar_product(vector(p2, p1), n2),
                            scalar_product(vector(p2, q1), n2),
                            scalar_product(vector(p2, r1), n2) };

  // All the vertices of tr1 must be on the same side of tr2
  // Coplanarity is tolerated, so '1' and '0' should be allowed, but not '1' and '-1'
  if(CGAL::sign(sps[0]) == - CGAL::sign(sps[1]) || CGAL::sign(sps[1]) == - CGAL::sign(sps[2]))
    return std::make_pair(0, false);

  std::for_each(sps.begin(), sps.end(), [](FT& v) { v = abs(v); });
  const auto min_pos = std::min_element(sps.begin(), sps.end());
  const std::size_t min_id = static_cast<std::size_t>(std::distance(sps.begin(), min_pos));

  if(sps[min_id] > 0)
    are_triangles_known_to_be_disjoint = true;

  const Point_3& x1 = vertex(tr1, static_cast<int>(min_id));

  if(CGAL::internal::on_left_of_triangle_edge(x1, n2, p2, q2, k) &&
     CGAL::internal::on_left_of_triangle_edge(x1, n2, q2, r2, k) &&
     CGAL::internal::on_left_of_triangle_edge(x1, n2, r2, p2, k))
  {
    // the projection of `x1` is inside the triangle
    return std::make_pair(CGAL::internal::squared_distance_to_plane(n2, vector(p2, x1), k), true);
  }

  return std::make_pair(0, false);
}

} // namespace internal
} // namespace Distance_3

namespace internal {

template <typename K>
typename K::FT
squared_distance(const typename K::Triangle_3& tr1,
                 const typename K::Triangle_3& tr2,
                 const K& k)
{
  typedef typename K::FT                                                  FT;

  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();

  // ideally just limits<FT>::infinity|max(), but it is not available for exact NTs...
  FT global_min_sqd = squared_distance(vertex(tr1, 0), vertex(tr2, 0));

  bool are_triangles_known_to_be_disjoint = false;
  std::pair<Distance_3::internal::Segment_3_Segment_3_Result<K>, bool> ss_res;
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j)
    {
      ss_res = Distance_3::internal::test_edge_pair(
                 vertex(tr1, i%3), vertex(tr1, (i+1)%3), vertex(tr1, (i+2)%3),
                 vertex(tr2, j%3), vertex(tr2, (j+1)%3), vertex(tr2, (j+2)%3), k,
                 global_min_sqd, are_triangles_known_to_be_disjoint);

      if(ss_res.second)
        return ss_res.first.squared_distance;
    }
  }

  // Failed to find a minimum between segment pairs, explore vertex-triangle distances

#if 1
  std::pair<FT, bool> pt_res =
      Distance_3::internal::test_vertex_triangle(tr1, tr2, k, are_triangles_known_to_be_disjoint);
  if(pt_res.second)
    return pt_res.first;

  pt_res = Distance_3::internal::test_vertex_triangle(tr2, tr1, k, are_triangles_known_to_be_disjoint);
  if(pt_res.second)
    return pt_res.first;

  if(are_triangles_known_to_be_disjoint)
    return global_min_sqd;
  else
    return 0;
#else // A tiny bit less efficient, but a lot clearer!
  // @todo does not handle degenerate inputs
  if(!are_triangles_known_to_be_disjoint && CGAL::do_intersect(tr1, tr2))
    return 0;

  FT sqd_p1 = CGAL::squared_distance(vertex(tr1, 0), tr2);
  FT sqd_q1 = CGAL::squared_distance(vertex(tr1, 1), tr2);
  FT sqd_r1 = CGAL::squared_distance(vertex(tr1, 2), tr2);
  FT sqd_p2 = CGAL::squared_distance(vertex(tr2, 0), tr1);
  FT sqd_q2 = CGAL::squared_distance(vertex(tr2, 1), tr1);
  FT sqd_r2 = CGAL::squared_distance(vertex(tr2, 2), tr1);

  const FT m = std::min({sqd_p1, sqd_q1, sqd_r1, sqd_p2, sqd_q2, sqd_r2});

  return m;
#endif
}

template <typename K>
typename K::Comparison_result
compare_squared_distance_disjoint(const typename K::Triangle_3& tr1,
                         const typename K::Triangle_3& tr2,
                         const typename K::FT squared_distance,
                         const K& k)
{
  typedef CGAL::Uncertain<Comparison_result> U;

  assert(!CGAL::do_intersect(tr1,tr2));

  auto compare_squared_distance=[k](const auto A, auto B, const typename K::FT squared_distance){
    return compare(internal::squared_distance(A,B,k), squared_distance);
  };

  auto triangle_contains_projected_point=[](const typename K::Triangle_3 &tr, const typename K::Point_3 &p){
    typename K::Point_3 p_project = tr.supporting_plane().projection(p);
    return tr.has_on(p_project);
  };

  U tri_tri_result(CGAL::LARGER);
  auto update_result=[]( U &res, const U u){
    res = U( min(res.inf(),u.inf()), min(res.sup(),u.sup()));
  };
  //First, we check the distance between each pair of segments
	//After, we check the distance between each vertex to the opposite triangle
  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j)
      //If the distance between the lines is more than 3, the distance between the segments is more than 3 and thus, we not compute it.
      //If the distance between the lines is less than 3, it is computed again in the distance between the segments but this gains a bit of time because
      //most pair of lines are at distance more than 3
      if(compare_squared_distance(typename K::Line_3(tr1[i],tr1[(i+1)%3]), typename K::Line_3(tr2[j],tr2[(j+1)%3]),squared_distance)!=CGAL::LARGER){
        U seg_seg_result(compare_squared_distance(typename K::Segment_3(tr1[i],tr1[(i+1)%3]), typename K::Segment_3(tr2[j],tr2[(j+1)%3]),squared_distance));
        if(seg_seg_result.is_certain() && (seg_seg_result.make_certain()==CGAL::SMALLER))
          return CGAL::SMALLER;
        update_result(tri_tri_result, seg_seg_result);
      }
    //If the distance to the plane is more than 3, the distance to the triangle is more than 3 and we exit.
    //We already test that all edges are at distance more than 3 so we can simply test if the projection of the vertex is inside the triangle.
    U pts_plane_result(compare_squared_distance(tr1.supporting_plane(),tr2[i],squared_distance));
    if(CGAL::possibly(pts_plane_result!=CGAL::LARGER) && triangle_contains_projected_point(tr1,tr2[i])){
      if(pts_plane_result.is_certain() && (pts_plane_result.make_certain()==CGAL::SMALLER))
          return CGAL::SMALLER;
        update_result(tri_tri_result, pts_plane_result);
    }
    pts_plane_result = compare_squared_distance(tr1[i],tr2.supporting_plane(),squared_distance);
    if(CGAL::possibly(pts_plane_result!=CGAL::LARGER) && triangle_contains_projected_point(tr2,tr1[i])){
      if(pts_plane_result.is_certain() && (pts_plane_result.make_certain()==CGAL::SMALLER))
          return CGAL::SMALLER;
        update_result(tri_tri_result, pts_plane_result);
    }
  }
  return tri_tri_result;
}

template <typename K>
typename K::Comparison_result
compare_squared_distance(const typename K::Triangle_3& tr1,
                         const typename K::Triangle_3& tr2,
                         const typename K::FT squared_distance,
                         const K& k)
{
  if(CGAL::do_intersect(tr1,tr2))
    return CGAL::is_zero(squared_distance)?CGAL::EQUAL:CGAL::SMALLER;
  return compare_squared_distance_disjoint(tr1,tr2,squared_distance,k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Triangle_3<K>& tr1,
                 const Triangle_3<K>& tr2)
{
  return K().compute_squared_distance_3_object()(tr1, tr2);
}


} // namespace CGAL

#endif // CGAL_DISTANCE_3_TRIANGLE_3_TRIANGLE_3_H
