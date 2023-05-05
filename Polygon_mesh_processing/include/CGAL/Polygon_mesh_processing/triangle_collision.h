// Copyright (c) 2023 INRIA Lorraine LORIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Valque Léo

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_COLLISION_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_COLLISION_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

namespace CGAL{
namespace Polygon_mesh_processing{

enum class MotionStep{X_MOVING, Y_MOVING, Z_MOVING, COPLANAR};
enum class CollisionType{EDGE_EDGE, VERTEX_TRIANGLE, VERTEX_VERTEX};

struct VertexVertexInfo{
	VertexVertexInfo(size_t p1_, size_t p2_):p1(p1_),p2(p2_){}
	size_t p1;
	size_t p2;
};

struct EdgeEdgeInfo{
	EdgeEdgeInfo(std::pair<size_t, size_t> ei1_, std::pair<size_t, size_t> ei2_): ei1(ei1_), ei2(ei2_){}
	std::pair<size_t, size_t> ei1;
	std::pair<size_t, size_t> ei2;
};

struct VertexTriangleInfo{
	VertexTriangleInfo(size_t pi_, size_t fi_): pi(pi_), fi(fi_){}
	size_t pi;
	size_t fi;
};

struct Collision{
	template<class Info>
	Collision(CollisionType type_, MotionStep during_, Info info_):type(type_), collideDuring(during_), info(info_){}

	CollisionType type;
	MotionStep collideDuring;
	std::variant<VertexVertexInfo, EdgeEdgeInfo, VertexTriangleInfo> info;
};

namespace internal{

template<int d>
struct Index{
	Index():elt(-1){}
	explicit Index(size_t elt_):elt(elt_){}
	operator size_t() const {return elt;}
	bool operator==(Index<d> other) const {return elt==other.elt;}
	size_t elt;
};

typedef Index<0>                                          P_index; //index of points vector
typedef Index<1>                                          F_index; //index of faces vector

struct Parameter_Face_Vertex{
	Parameter_Face_Vertex(F_index input) : is_face_parameter(true), face_parameter(input){}
	Parameter_Face_Vertex(P_index input) : is_face_parameter(false), vertex_parameter(input){}

	bool is_face_parameter;
	F_index face_parameter;
	P_index vertex_parameter;
};

typedef std::pair< P_index, P_index > Parameter_Edge;

typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3, std::vector< Parameter_Face_Vertex >::iterator > Face_Vertex_Box;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3, std::set< Parameter_Edge>::iterator> Edge_Box;

template<class K>
struct Functor_points_collide{
	Functor_points_collide(std::vector< typename K::Point_3> *points_, 
	                       std::vector< std::vector<size_t>> *triangles_, 
						   std::vector< Collision > *out_collision_,
	                       std::vector< typename K::Point_3 > *points_move_x_,
						   std::vector< typename K::Point_3 > *points_move_y_,
						   std::vector< typename K::Point_3 > *points_move_z_)
						  :points(points_),
						   triangles(triangles_), 
						   out_collision(out_collision_),
						   points_move_x(points_move_x_),
						   points_move_y(points_move_y_),
						   points_move_z(points_move_z_){}
	
	void operator()(const Face_Vertex_Box& a_box, const Face_Vertex_Box& b_box) {		
		P_index pa = a_box.handle()->vertex_parameter;
		P_index pb = b_box.handle()->vertex_parameter;

		MotionStep collideDuring;
		if( does_vertices_collide(pa, pb, collideDuring)){
			out_collision->emplace_back(CollisionType::VERTEX_VERTEX, collideDuring, VertexVertexInfo(pa,pb));
		}
	}

	bool does_vertices_collide(size_t pia, size_t pib, MotionStep &collideDuring);

	std::vector< typename K::Point_3 > *points_move_x;
	std::vector< typename K::Point_3 > *points_move_y;
	std::vector< typename K::Point_3 > *points_move_z;
	std::vector< Collision > *out_collision;
	std::vector< typename K::Point_3 > *points;
	std::vector< std::vector<size_t>> *triangles;
};

template<class K>
struct Functor_edge_collide{
	Functor_edge_collide(std::vector< typename K::Point_3> *points_, 
	                       std::vector< std::vector<size_t>> *triangles_, 
						   std::vector< Collision > *out_collision_,
	                       std::vector< typename K::Point_3 > *points_move_x_,
						   std::vector< typename K::Point_3 > *points_move_y_,
						   std::vector< typename K::Point_3 > *points_move_z_, 
						   bool simplicialComplexOption_)
						  :points(points_),
						   triangles(triangles_), 
						   out_collision(out_collision_),
						   points_move_x(points_move_x_),
						   points_move_y(points_move_y_),
						   points_move_z(points_move_z_), 
						   simplicialComplexOption(simplicialComplexOption_){}
	
	void operator()(const Edge_Box& a_box, const Edge_Box& b_box) {		
		Parameter_Edge a = *(a_box.handle());
		Parameter_Edge b = *(b_box.handle());
		
		MotionStep collideDuring;
		if(does_edge_edge_collide(a.first, a.second, b.first, b.second, collideDuring, simplicialComplexOption)){
			out_collision->emplace_back(CollisionType::EDGE_EDGE, collideDuring, EdgeEdgeInfo(std::pair<size_t,size_t> (a.first, a.second),
			                                                                                          std::pair<size_t,size_t> (b.first, b.second)));
		}
	}

	bool does_edge_edge_collide(size_t a1, size_t a2, size_t b1, size_t b2, MotionStep &collideDuring, bool simplicialComplexOption);

	
	std::vector< typename K::Point_3 > *points_move_x;
	std::vector< typename K::Point_3 > *points_move_y;
	std::vector< typename K::Point_3 > *points_move_z;
	std::vector< Collision > *out_collision;
	std::vector< typename K::Point_3 > *points;
	std::vector< std::vector<size_t>> *triangles;
	bool simplicialComplexOption;
};

template<class K>
struct Functor_face_collide{
	Functor_face_collide(std::vector< typename K::Point_3> *points_, 
	                       std::vector< std::vector<size_t>> *triangles_, 
						   std::vector< Collision > *out_collision_,
	                       std::vector< typename K::Point_3 > *points_move_x_,
						   std::vector< typename K::Point_3 > *points_move_y_,
						   std::vector< typename K::Point_3 > *points_move_z_, 
						   bool simplicialComplexOption_)
						  :points(points_),
						   triangles(triangles_), 
						   out_collision(out_collision_),
						   points_move_x(points_move_x_),
						   points_move_y(points_move_y_),
						   points_move_z(points_move_z_), 
						   simplicialComplexOption(simplicialComplexOption_){}

	void operator()(const Face_Vertex_Box& p_box, const Face_Vertex_Box& tr_box){
		std::vector<size_t> &tr = (*triangles)[tr_box.handle()->face_parameter];
		size_t pi = p_box.handle()->vertex_parameter;

		MotionStep collideDuring;
		if( does_vertex_triangle_collide(tr[0],tr[1], tr[2], pi, collideDuring, simplicialComplexOption)){
			out_collision->emplace_back(CollisionType::VERTEX_TRIANGLE, collideDuring, VertexTriangleInfo(pi, tr_box.handle()->face_parameter));
		}
	}

	bool does_vertex_triangle_collide(size_t t1, size_t t2, size_t t3, size_t v, MotionStep &collideDuring, bool simplicialComplexOption);
	
	std::vector< typename K::Point_3 > *points_move_x;
	std::vector< typename K::Point_3 > *points_move_y;
	std::vector< typename K::Point_3 > *points_move_z;
	std::vector< Collision > *out_collision;
	std::vector< typename K::Point_3 > *points;
	std::vector< std::vector<size_t>> *triangles;
	bool simplicialComplexOption;
};

template<class K>
bool seg_do_intersect(const typename K::Segment_2 s1, const typename K::Segment_2 s2, bool simplicialComplexOption){
	if(simplicialComplexOption){
		//if(CGAL::do_intersect(s1, s2.start()) || CGAL::do_intersect(s1, s2.end()) || CGAL::do_intersect(s2, s1.start()) || CGAL::do_intersect(s2, s1.end()) || CGAL::do_intersect(s1,s2))
		return CGAL::do_intersect(s1,s2);
		 
	} else if(CGAL::do_intersect(s1, s2.start()) || CGAL::do_intersect(s1, s2.end()) || CGAL::do_intersect(s2, s1.start()) || CGAL::do_intersect(s2, s1.end()) ){
			return false;
		} else {
			return CGAL::do_intersect(s1,s2);
		}
		
}

template<class K>
bool seg_do_intersect(const typename K::Segment_3 s1, const typename K::Segment_3 s2, bool simplicialComplexOption){
	
    if(simplicialComplexOption){
		if(s1.is_degenerate())
			if(s2.is_degenerate())
				return false;
			else
				return CGAL::do_intersect(s2, s1.start());
		else if(s2.is_degenerate())
			return CGAL::do_intersect(s1, s2.start());
		else
			return CGAL::do_intersect(s1,s2);
	} else {
		if(s1.is_degenerate() || s2.is_degenerate())
			return false;
		else if(CGAL::do_intersect(s1, s2.start()) || CGAL::do_intersect(s1, s2.end()) || CGAL::do_intersect(s2, s1.start()) || CGAL::do_intersect(s2, s1.end()) )
			return false;
		else
			return CGAL::do_intersect(s1,s2);
	}
}

template<class K>
bool Functor_points_collide<K>::does_vertices_collide(size_t pia, size_t pib, MotionStep &collideDuring){
  auto on_floor=[](typename K::Point_3 &p){
	return typename K::Point_2(p.x(), p.y());
  };
  auto on_backwall=[](typename K::Point_3 &p){
	return typename K::Point_2(p.x(), p.z());
  };	
  auto on_sidewall=[](typename K::Point_3 &p){
	return typename K::Point_2(p.y(), p.z());
  };	

  std::array<typename K::Point_3, 4> pa = {(*points)[pia],(*points_move_x)[pia],(*points_move_y)[pia],(*points_move_z)[pia]};
  std::array<typename K::Point_3, 4> pb = {(*points)[pib],(*points_move_x)[pib],(*points_move_y)[pib],(*points_move_z)[pib]};
  typename K::Equal_2 equal;

  //move on x
  if(equal(on_sidewall(pa[0]),on_sidewall(pb[0]))){
	auto comp_bef=compare_x(pa[0],pb[0]);
	auto comp_aft=compare_x(pa[1],pb[1]);
	if(comp_bef!=comp_aft && comp_aft!=EQUAL){
		collideDuring=MotionStep::X_MOVING;
  		return true;
	}  
  }
  
  //move on y
  if(equal(on_backwall(pa[1]),on_backwall(pb[1]))){
	auto comp_bef=compare_y(pa[1],pb[1]);
	auto comp_aft=compare_y(pa[2],pb[2]);
	if(comp_bef!=comp_aft && comp_aft!=EQUAL){
		collideDuring=MotionStep::Y_MOVING;
  		return true;
	}  
  }
  
  //move on z
  if(equal(on_floor(pa[2]),on_floor(pb[2]))){
	auto comp_bef=compare_z(pa[2],pb[2]);
	auto comp_aft=compare_z(pa[3],pb[3]);
	if(comp_bef!=comp_aft && comp_aft!=EQUAL){
		collideDuring=MotionStep::Z_MOVING;
  		return true;
	}  
  }
  return false;
}

template<class K>
bool Functor_face_collide<K>::does_vertex_triangle_collide(size_t t1, size_t t2, size_t t3, size_t v, MotionStep &collideDuring, bool simplicialComplexOption){

  if(t1==t2 || t1==t3 || t2==t3 || t1==v || t2==v || t3==v)
  	return false;

  std::array< typename K::Point_3, 4> p = {(*points)[t1],(*points)[t2],(*points)[t3],(*points)[v]};
  std::array< typename K::Point_3, 4> px = {(*points_move_x)[t1],(*points_move_x)[t2],(*points_move_x)[t3],(*points_move_x)[v]};
  std::array< typename K::Point_3, 4> py = {(*points_move_y)[t1],(*points_move_y)[t2],(*points_move_y)[t3],(*points_move_y)[v]};
  std::array< typename K::Point_3, 4> pz = {(*points_move_z)[t1],(*points_move_z)[t2],(*points_move_z)[t3],(*points_move_z)[v]};

  CGAL::Orientation orientation_begin =   CGAL::orientation(p[0],p[1],p[2],p[3]);
  CGAL::Orientation orientation_x = CGAL::orientation(px[0],px[1],px[2],px[3]);
  CGAL::Orientation orientation_y = CGAL::orientation(py[0],py[1],py[2],py[3]);
  CGAL::Orientation orientation_z = CGAL::orientation(pz[0],pz[1],pz[2],pz[3]);
  
  typename K::Triangle_3 tr(p[0], p[1], p[2]);
  typename K::Triangle_2 tr_on_x (typename K::Point_2(p[0].y(),p[0].z())  ,typename K::Point_2(p[1].y(),p[1].z())  ,typename K::Point_2(p[2].y(),p[2].z()));
  typename K::Triangle_2 trx_on_y(typename K::Point_2(px[0].x(),px[0].z()),typename K::Point_2(px[1].x(),px[1].z()),typename K::Point_2(px[2].x(),px[2].z()));
  typename K::Triangle_2 try_on_z(typename K::Point_2(py[0].x(),py[0].y()),typename K::Point_2(py[1].x(),py[1].y()),typename K::Point_2(py[2].x(),py[2].y()));

  //move on x
  if( orientation_begin != COPLANAR &&
	  orientation_x != COPLANAR &&
  	  orientation_x != orientation_begin &&
  	  tr_on_x.bounded_side(typename K::Point_2(p[3].y(),p[3].z())) != ON_UNBOUNDED_SIDE){
	  collideDuring=MotionStep::X_MOVING;
	  return true;
  }
  
  //move on y
  if(orientation_y != COPLANAR &&
     orientation_y != orientation_x &&
     trx_on_y.bounded_side( typename K::Point_2(px[3].x(),px[3].z())) != ON_UNBOUNDED_SIDE &&
     (orientation_x != COPLANAR || 
	 (orientation_begin != COPLANAR && orientation_begin != orientation_y))){
	collideDuring=MotionStep::Y_MOVING;
	return true;
  }
  
  //move on z
  if(orientation_z != CGAL::COPLANAR &&
  	orientation_z != orientation_y &&
  	try_on_z.bounded_side(typename K::Point_2(py[3].x(),py[3].y())) != ON_UNBOUNDED_SIDE &&
  	(orientation_y != CGAL::COPLANAR || 
		(orientation_x != CGAL::COPLANAR && orientation_x != orientation_z) || 
		(orientation_x == CGAL::COPLANAR && orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_z))){
	collideDuring=MotionStep::Z_MOVING;
  	return true;
  }	
  	
  //move epsilon
  if(simplicialComplexOption &&
    orientation_z == COPLANAR &&
  	typename K::Triangle_3(pz[0], pz[1], pz[2]).has_on(pz[3]) &&
  	!( pz[3] == pz[2] || pz[3] == pz[1] || pz[3] == pz[0])){
	collideDuring=MotionStep::COPLANAR;
  	return true;
  }//*/
  
  return false;
}

template<class K>
bool Functor_edge_collide<K>::does_edge_edge_collide(size_t a1, size_t a2, size_t b1, size_t b2, MotionStep &collideDuring, bool simplicialComplexOption){
  if(a1==b1 || a1==b2 || a2==b1 || a2==b2)
  	return false; //they can't intersect if they are a common vertex

  std::array<typename K::Point_3, 4> p = {(*points)[a1],(*points)[a2],(*points)[b1],(*points)[b2]};
  std::array<typename K::Point_3, 4> px = {(*points_move_x)[a1],(*points_move_x)[a2],(*points_move_x)[b1],(*points_move_x)[b2]};
  std::array<typename K::Point_3, 4> py = {(*points_move_y)[a1],(*points_move_y)[a2],(*points_move_y)[b1],(*points_move_y)[b2]};
  std::array<typename K::Point_3, 4> pz = {(*points_move_z)[a1],(*points_move_z)[a2],(*points_move_z)[b1],(*points_move_z)[b2]};

  CGAL::Orientation orientation_begin   = orientation(p[0],p[1],p[2],p[3]);
  CGAL::Orientation orientation_x = orientation(px[0],px[1],px[2],px[3]);
  CGAL::Orientation orientation_y = orientation(py[0],py[1],py[2],py[3]);
  CGAL::Orientation orientation_z = orientation(pz[0],pz[1],pz[2],pz[3]);

  //move on x
  if( orientation_begin != CGAL::COPLANAR &&
  	  orientation_x != CGAL::COPLANAR &&
  	  orientation_x != orientation_begin)
  {
	typename K::Segment_2 a_on_x(typename K::Point_2(p[0].y(),p[0].z()), typename K::Point_2(p[1].y(), p[1].z()));
  	typename K::Segment_2 b_on_x(typename K::Point_2(p[2].y(),p[2].z()), typename K::Point_2(p[3].y(), p[3].z()));
	if(seg_do_intersect<K>(a_on_x,b_on_x,simplicialComplexOption))
	{
		collideDuring=MotionStep::X_MOVING;
		return true;  
	}
  }

  //move on y
  if(orientation_y != CGAL::COPLANAR &&
     orientation_y != orientation_x &&
     (orientation_x != CGAL::COPLANAR || (orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_y)))
  {
    typename K::Segment_2 ax_on_y(typename K::Point_2(px[0].x(),px[0].z()), typename K::Point_2(px[1].x(), px[1].z()));
    typename K::Segment_2 bx_on_y(typename K::Point_2(px[2].x(),px[2].z()), typename K::Point_2(px[3].x(), px[3].z()));
	if(seg_do_intersect<K>(ax_on_y,bx_on_y,simplicialComplexOption))
	{
		collideDuring=MotionStep::Y_MOVING;
		return true;
	}
  }

  //move on z
  if(orientation_z != CGAL::COPLANAR &&
  	 orientation_z != orientation_y &&
  	(orientation_y != CGAL::COPLANAR || 
		(orientation_x != CGAL::COPLANAR && orientation_x != orientation_z) || 
		(orientation_x == CGAL::COPLANAR && orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_z)))
  {

	typename K::Segment_2 ay_on_z(typename K::Point_2(py[0].x(),py[0].y()), typename K::Point_2(py[1].x(), py[1].y()));
	typename K::Segment_2 by_on_z(typename K::Point_2(py[2].x(),py[2].y()), typename K::Point_2(py[3].x(), py[3].y()));
	if(seg_do_intersect<K>(ay_on_z,by_on_z,simplicialComplexOption))
	{
		collideDuring=MotionStep::Z_MOVING;
		return true;
	}
  }

  //move epsilon
  if(orientation_z == CGAL::COPLANAR &&
  	!( pz[2] == pz[0] || pz[2] == pz[1] || pz[3] == pz[1] || pz[3] == pz[0]))
  {
	if(seg_do_intersect<K>(typename K::Segment_3(pz[0],pz[1]), typename K::Segment_3(pz[2],pz[3]),simplicialComplexOption))
	{
		typename K::Segment_2 ay_on_z(typename K::Point_2(py[0].x(),py[0].y()), typename K::Point_2(py[1].x(), py[1].y()));
		typename K::Segment_2 by_on_z(typename K::Point_2(py[2].x(),py[2].y()), typename K::Point_2(py[3].x(), py[3].y()));
		collideDuring=MotionStep::COPLANAR;
		return true;
	}
  }

  return false;
}


}

//Explanation of the moving sheme:
//
//Disturb symbolically the points using SOS, we can assume that the two segments are not coplanar during motion except possibly at the end
//this disturbtion and the epsilon is not in the code, it's just in the proof to justify that we can go through the coplanar case except at the end
//move the full motion minus epsilon on x and check collision
//move the full motion minus epsilon on y and check collision
//move the full motion minus epsilon on z and check collision
//move the remaining motion, concretely it corresponds to check the coplanar case at the end

template<class K,
         class Motion>
bool detect_collisions_during_motion(std::vector<typename K::Point_3> &points, std::vector< std::vector<size_t> > &triangles, const Motion &motion, std::vector<Collision> &out_collision, bool simplicialComplexOption){
	auto motion_x=[motion](const typename K::Point_3 &p, size_t pi){
		return typename K::Point_3(motion(p.x(),pi),p.y(),p.z());
	};
	auto motion_xy=[motion](const typename K::Point_3 &p, size_t pi){
		return typename K::Point_3(motion(p.x(),pi),motion(p.y(),pi),p.z());
	};
	auto motion_xyz=[motion](const typename K::Point_3 &p, size_t pi){
		return typename K::Point_3(motion(p.x(),pi),motion(p.y(),pi),motion(p.z(),pi));
	};

	auto extending_bbox=[motion, motion_xyz](typename K::Point_3 &p, size_t pi){
		return p.bbox()+motion_xyz(p, pi).bbox();
	};

	auto triangle_bbox=[motion, extending_bbox](std::vector<size_t> &tr, std::vector<typename K::Point_3> &points){
		return extending_bbox(points[tr[0]],tr[0])+extending_bbox(points[tr[1]],tr[1])+extending_bbox(points[tr[2]],tr[2]);
	};
	
	std::vector< internal::Parameter_Face_Vertex > face_parameters;
	std::set< internal::Parameter_Edge > edge_parameters;
	std::vector< internal::Parameter_Face_Vertex > vertex_parameters;

	std::vector< typename K::Point_3 > points_move_x;
	std::vector< typename K::Point_3 > points_move_y;
	std::vector< typename K::Point_3 > points_move_z;
	for(size_t i=0; i<points.size(); ++i){
		auto &p = points[i];
		points_move_x.push_back(motion_x(p, i));
		points_move_y.push_back(motion_xy(p, i));
		points_move_z.push_back(motion_xyz(p, i));
	}

	for(size_t i=0; i<points.size(); ++i){
		vertex_parameters.emplace_back(internal::P_index(i));
	}
	for(size_t i=0; i<triangles.size(); ++i){
		face_parameters.emplace_back(internal::F_index(i));
		assert(triangles[i].size()==3);
	}

	for(auto &t : triangles){
		if(t[0]<t[1])
			edge_parameters.emplace(internal::P_index(t[0]), internal::P_index(t[1]));
		else
			edge_parameters.emplace(internal::P_index(t[1]), internal::P_index(t[0]));
		if(t[0]<t[2])
			edge_parameters.emplace(internal::P_index(t[0]), internal::P_index(t[2]));
		else
			edge_parameters.emplace(internal::P_index(t[2]), internal::P_index(t[0]));
		if(t[1]<t[2])
			edge_parameters.emplace(internal::P_index(t[1]), internal::P_index(t[2]));
		else
			edge_parameters.emplace(internal::P_index(t[2]), internal::P_index(t[1]));
	}

	std::vector< internal::Face_Vertex_Box > face_boxes;
	std::vector< internal::Edge_Box > edge_boxes;
	std::vector< internal::Face_Vertex_Box > vertex_boxes;
	
	for(auto it=face_parameters.begin(); it!=face_parameters.end(); ++it){
		face_boxes.emplace_back(triangle_bbox(triangles[it->is_face_parameter], points), it);
	}
	for(auto it=vertex_parameters.begin(); it!=vertex_parameters.end(); ++it)
		vertex_boxes.emplace_back(extending_bbox(points[it->vertex_parameter], it->vertex_parameter), it);
	
	for(auto it=edge_parameters.begin(); it!=edge_parameters.end(); ++it)
		edge_boxes.emplace_back(extending_bbox(points[it->first], it->first)+extending_bbox(points[it->second], it->second), it);
		
	internal::Functor_edge_collide<K> functor1(&points, &triangles, &out_collision, &points_move_x, &points_move_y, &points_move_z, simplicialComplexOption);
    CGAL::box_self_intersection_d( edge_boxes.begin(), edge_boxes.end(), functor1);
    internal::Functor_face_collide<K> functor2(&points, &triangles, &out_collision,  &points_move_x, &points_move_y, &points_move_z, simplicialComplexOption);
    CGAL::box_intersection_d( vertex_boxes.begin(), vertex_boxes.end(), face_boxes.begin(), face_boxes.end(), functor2);
	internal::Functor_points_collide<K> functor3(&points, &triangles, &out_collision,  &points_move_x, &points_move_y, &points_move_z);
    CGAL::box_self_intersection_d( vertex_boxes.begin(), vertex_boxes.end(),functor3);

    return out_collision.size()!=0;
}

}
}

#endif