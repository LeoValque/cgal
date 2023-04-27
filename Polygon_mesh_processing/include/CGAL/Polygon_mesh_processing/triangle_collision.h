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

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef CGAL::Bbox_2       Bbox_2;
typedef CGAL::Bbox_3       Bbox_3;

typedef Kernel::Vector_2   Vector_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Segment_2  Segment_2;
typedef Kernel::Line_2     Line_2;
typedef Kernel::Point_2    Point_2;

typedef Kernel::Point_3    Point_3;

typedef Kernel::Vector_3   Vector_3;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Segment_3  Segment_3;
typedef Kernel::Line_3     Line_3;
typedef Kernel::Plane_3    Plane_3;

enum class MotionStep{X_MOVING, Y_MOVING, Z_MOVING, COPLANAR};

struct Intersection{
	Intersection(size_t pi_, size_t fi_, MotionStep intersectDuring_): pi(pi_), fi(fi_), intersectionDuring(intersectDuring_), is_edge_edge(false){}
	Intersection(std::pair<size_t, size_t> ei1_, std::pair<size_t, size_t> ei2_, MotionStep intersectDuring_): ei1(ei1_), ei2(ei2_), intersectionDuring(intersectDuring_), is_edge_edge(true){}
	//Intersection(size_t p1, size_t p2, size_t p3, size_t p4, MotionStep intersectDuring_, bool is_edge_edge_):feature_points{p1,p2,p3,p4},is_edge_edge(is_edge_edge_),intersectionDuring(intersectDuring_){}

	bool is_edge_edge;
	size_t pi;
	size_t fi;
	std::pair<size_t, size_t> ei1;
	std::pair<size_t, size_t> ei2;
	//std::array<size_t, 4> feature_points;
	MotionStep intersectionDuring;

};

struct FunctorSnap{
	virtual Kernel::FT snap_FT(Kernel::FT x, size_t pi)=0;
	Point_3 snap_x(const Point_3 &p, size_t pi){
		return Point_3(snap_FT(p.x(),pi),p.y(),p.z());
	}
	Point_3 snap_xy(const Point_3 &p, size_t pi){
		return Point_3(snap_FT(p.x(),pi),snap_FT(p.y(),pi),p.z());
	}
	Point_3 snap_xyz(const Point_3 &p, size_t pi){
		return Point_3(snap_FT(p.x(),pi),snap_FT(p.y(),pi),snap_FT(p.z(),pi));
	}
};

Bbox_3 extending_bbox(Point_3 &p, FunctorSnap &snap);


bool does_triangles_collide(std::vector<Point_3> &points, std::vector< std::vector<size_t> > &triangles, FunctorSnap &snap, std::vector<Intersection> &out_intersection, bool SimplicialComplexOption=false);

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/box_intersection_d.h>
#include <filesystem>

#define VERBOSE_OUT if(verbose) std::cout << "Debug: " << std::filesystem::path(__FILE__).filename().string() << ":" << __LINE__ << " " << __func__ <<std::endl

typedef std::chrono::_V2::high_resolution_clock high_resolution_clock;
typedef std::chrono::microseconds microseconds;
typedef std::chrono::seconds seconds;

const bool verbose=false;
//const bool option_simplicial_complex=true;
//const bool optionFaceNull=true;

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

size_t num_pair=0;

Kernel::FT two(2);

Bbox_3 extending_bbox(Point_3 &p, size_t pi, FunctorSnap &snap){
	return p.bbox()+snap.snap_xyz(p, pi).bbox();
}

Bbox_3 bbox(std::vector<size_t> &tr, std::vector<Point_3> &points, FunctorSnap &snap){
	return extending_bbox(points[tr[0]],tr[0],snap)+extending_bbox(points[tr[1]],tr[1],snap)+extending_bbox(points[tr[2]],tr[2],snap);
}

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

struct Functor_edge_intersect{
	Functor_edge_intersect(std::vector<Point_3> *points_, std::vector< std::vector<size_t>> *triangles_, std::vector<Intersection> *out_intersection_,
	                      std::vector< Point_3 > *points_round_x_,std::vector< Point_3 > *points_round_y_,std::vector< Point_3 > *points_round_z_, bool simplicialComplexOption_)
						  :points(points_),triangles(triangles_), out_intersection(out_intersection_),
						   points_round_x(points_round_x_),points_round_y(points_round_y_),points_round_z(points_round_z_), simplicialComplexOption(simplicialComplexOption_){}
	
	void operator()(const Edge_Box& a_box, const Edge_Box& b_box) {
		//ListIntersection &list_intersections = *list_intersections_;
		
		Parameter_Edge a = *(a_box.handle());
		Parameter_Edge b = *(b_box.handle());
		
		bool is_coplanar=false;
		MotionStep intersectDuring;
		if(edges_do_intersect(a.first, a.second, b.first, b.second, intersectDuring, simplicialComplexOption)){
			out_intersection->emplace_back(std::pair<size_t,size_t> (a.first, a.second), std::pair<size_t,size_t> (b.first, b.second), intersectDuring);
		}
	}

	bool edges_do_intersect(size_t a1, size_t a2, size_t b1, size_t b2, MotionStep &intersectDuring, bool simplicialComplexOption);

	
	std::vector< Point_3 > *points_round_x;
	std::vector< Point_3 > *points_round_y;
	std::vector< Point_3 > *points_round_z;
	std::vector<Intersection> *out_intersection;
	std::vector<Point_3> *points;
	std::vector< std::vector<size_t>> *triangles;
	bool simplicialComplexOption;
};

struct Functor_face_intersect{
	Functor_face_intersect(std::vector<Point_3> *points_, std::vector< std::vector<size_t>> *triangles_, std::vector<Intersection> *out_intersection_,
	                       std::vector< Point_3 > *points_round_x_,std::vector< Point_3 > *points_round_y_,std::vector< Point_3 > *points_round_z_, bool simplicialComplexOption_)
						   :points(points_),triangles(triangles_), out_intersection(out_intersection_),
						points_round_x(points_round_x_),points_round_y(points_round_y_),points_round_z(points_round_z_), simplicialComplexOption(simplicialComplexOption_){}

	void operator()(const Face_Vertex_Box& a_it, const Face_Vertex_Box& b_it){
		//ListIntersection &list_intersections = *list_intersections_;
		std::vector<size_t> &b = (*triangles)[b_it.handle()->face_parameter];
		size_t pi = a_it.handle()->vertex_parameter;

		MotionStep intersectDuring;
		if( faces_do_intersect(b[0],b[1], b[2], pi, intersectDuring, simplicialComplexOption)){
			//out_intersection->emplace_back(b[0], b[1], b[2], pi, intersectDuring,false);
			out_intersection->emplace_back(pi, b_it.handle()->face_parameter, intersectDuring);
		}
	}

	bool faces_do_intersect(size_t t1, size_t t2, size_t t3, size_t v, MotionStep &intersectDuring, bool simplicialComplexOption);
	
	std::vector< Point_3 > *points_round_x;
	std::vector< Point_3 > *points_round_y;
	std::vector< Point_3 > *points_round_z;
	std::vector<Intersection> *out_intersection;
	std::vector<Point_3> *points;
	std::vector< std::vector<size_t>> *triangles;
	bool simplicialComplexOption;
};

bool seg_do_intersect(const Segment_2 s1, const Segment_2 s2, bool simplicialComplexOption){
	if(simplicialComplexOption){
		//if(CGAL::do_intersect(s1, s2.start()) || CGAL::do_intersect(s1, s2.end()) || CGAL::do_intersect(s2, s1.start()) || CGAL::do_intersect(s2, s1.end()) || CGAL::do_intersect(s1,s2))
		return CGAL::do_intersect(s1,s2);
		 
	} else if(CGAL::do_intersect(s1, s2.start()) || CGAL::do_intersect(s1, s2.end()) || CGAL::do_intersect(s2, s1.start()) || CGAL::do_intersect(s2, s1.end()) ){
			return false;
		} else {
			return CGAL::do_intersect(s1,s2);
		}
		
}

bool seg_do_intersect(const Segment_3 s1, const Segment_3 s2, bool simplicialComplexOption){
	
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

/*
bool points_do_intersect(Point_3 a, Point_3 b, FunctorSnap &snap){
  std::array<Point_3, 4> pa = {a,snap.snap_x(a),snap.snap_xy(a),snap.snap_xyz(a)};
  std::array<Point_3, 4> pb = {b,snap.snap_x(b),snap.snap_xy(b),snap.snap_xyz(b)};
  Kernel::Equal_2 equal;

  //round on x
  if(equal(on_sidewall(pa[0]),on_sidewall(pb[0]))){
	auto comp_bef=CGAL::compare_x(pa[0],pb[0]);
	auto comp_aft=CGAL::compare_x(pa[1],pb[1]);
	if(comp_bef!=comp_aft && comp_aft!=CGAL::EQUAL){
		VERBOSE_OUT;
  		return true;
	}  
  }
  
  //round on y
  if(equal(on_backwall(pa[1]),on_backwall(pb[1]))){
	auto comp_bef=CGAL::compare_y(pa[1],pb[1]);
	auto comp_aft=CGAL::compare_y(pa[2],pb[2]);
	if(comp_bef!=comp_aft && comp_aft!=CGAL::EQUAL){
		VERBOSE_OUT;
  		return true;
	}  
  }
  
  //round on z
  if(equal(on_floor(pa[2]),on_floor(pb[2]))){
	auto comp_bef=CGAL::compare_z(pa[2],pb[2]);
	auto comp_aft=CGAL::compare_z(pa[3],pb[3]);
	if(comp_bef!=comp_aft && comp_aft!=CGAL::EQUAL){
		VERBOSE_OUT;
  		return true;
	}  
  }
  return false;
}*/

bool Functor_face_intersect::faces_do_intersect(size_t t1, size_t t2, size_t t3, size_t v, MotionStep &intersectDuring, bool simplicialComplexOption){

  if(t1==t2 || t1==t3 || t2==t3 || t1==v || t2==v || t3==v)
  	return false;

  std::array<Point_3, 4> p = {(*points)[t1],(*points)[t2],(*points)[t3],(*points)[v]};
  std::array<Point_3, 4> px = {(*points_round_x)[t1],(*points_round_x)[t2],(*points_round_x)[t3],(*points_round_x)[v]};
  std::array<Point_3, 4> py = {(*points_round_y)[t1],(*points_round_y)[t2],(*points_round_y)[t3],(*points_round_y)[v]};
  std::array<Point_3, 4> pz = {(*points_round_z)[t1],(*points_round_z)[t2],(*points_round_z)[t3],(*points_round_z)[v]};
  //std::array<Point_3, 4> px = {snap.snap_x(t1), snap.snap_x(t2), snap.snap_x(t3), snap.snap_x(v)};
  //std::array<Point_3, 4> py = {snap.snap_xy(t1), snap.snap_xy(t2), snap.snap_xy(t3), snap.snap_xy(v)};
  //std::array<Point_3, 4> pz = {snap.snap_xyz(t1), snap.snap_xyz(t2), snap.snap_xyz(t3), snap.snap_xyz(v)};

  CGAL::Orientation orientation_begin =   CGAL::orientation(p[0],p[1],p[2],p[3]);
  CGAL::Orientation orientation_round_x = CGAL::orientation(px[0],px[1],px[2],px[3]);
  CGAL::Orientation orientation_round_y = CGAL::orientation(py[0],py[1],py[2],py[3]);
  CGAL::Orientation orientation_round_z = CGAL::orientation(pz[0],pz[1],pz[2],pz[3]);
  
  Triangle_3 tr(p[0], p[1], p[2]);
  Triangle_2 tr_on_x(Point_2(p[0].y(),p[0].z()),Point_2(p[1].y(),p[1].z()),Point_2(p[2].y(),p[2].z()));
  Triangle_2 trx_on_y(Point_2(px[0].x(),px[0].z()),Point_2(px[1].x(),px[1].z()),Point_2(px[2].x(),px[2].z()));
  Triangle_2 try_on_z(Point_2(py[0].x(),py[0].y()),Point_2(py[1].x(),py[1].y()),Point_2(py[2].x(),py[2].y()));

  //round on x
  if( orientation_begin != CGAL::COPLANAR &&
	  orientation_round_x != CGAL::COPLANAR &&
  	  orientation_round_x != orientation_begin &&
  	  tr_on_x.bounded_side(Point_2(p[3].y(),p[3].z())) != CGAL::ON_UNBOUNDED_SIDE){
	  intersectDuring=MotionStep::X_MOVING;
	  return true;
  }
  
  //round on y
  if(orientation_round_y != CGAL::COPLANAR &&
     orientation_round_y != orientation_round_x &&
     trx_on_y.bounded_side(Point_2(px[3].x(),px[3].z())) != CGAL::ON_UNBOUNDED_SIDE &&
     (orientation_round_x != CGAL::COPLANAR || 
	 (orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_round_y))){
	intersectDuring=MotionStep::Y_MOVING;
	return true;
  }
  
  //round on z
  if(orientation_round_z != CGAL::COPLANAR &&
  	orientation_round_z != orientation_round_y &&
  	try_on_z.bounded_side(Point_2(py[3].x(),py[3].y())) != CGAL::ON_UNBOUNDED_SIDE &&
  	(orientation_round_y != CGAL::COPLANAR || 
		(orientation_round_x != CGAL::COPLANAR && orientation_round_x != orientation_round_z) || 
		(orientation_round_x == CGAL::COPLANAR && orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_round_z))){
	intersectDuring=MotionStep::Z_MOVING;
  	return true;
  }	
  	
  //round epsilon
  if(simplicialComplexOption &&
    orientation_round_z == CGAL::COPLANAR &&
  	Triangle_3(pz[0], pz[1], pz[2]).has_on(pz[3]) &&
  	!( pz[3] == pz[2] || pz[3] == pz[1] || pz[3] == pz[0])){
	intersectDuring=MotionStep::COPLANAR;
  	return true;
  }//*/
  
  return false;
}

bool Functor_edge_intersect::edges_do_intersect(size_t a1, size_t a2, size_t b1, size_t b2, MotionStep &intersectDuring, bool simplicialComplexOption){
  if(a1==b1 || a1==b2 || a2==b1 || a2==b2)
  	return false; //they can't intersect if they are a common vertex

  //TODO : il faut changer le egal avec l'optim sur l'arrondi sur float
  /*std::array<Point_3, 4> p = {a1,a2,b1,b2};
  std::array<Point_3, 4> px = {snap.snap_x(a1), snap.snap_x(a2), snap.snap_x(b1), snap.snap_x(b2)};
  std::array<Point_3, 4> py = {snap.snap_xy(a1), snap.snap_xy(a2), snap.snap_xy(b1), snap.snap_xy(b2)};
  std::array<Point_3, 4> pz = {snap.snap_xyz(a1), snap.snap_xyz(a2), snap.snap_xyz(b1), snap.snap_xyz(b2)};*/
  //num_pair++;
  //if((num_pair%10000)==0)
	//std::cout << num_pair << "\n";
  //std::cout << points->size() << " " << points_round_x->size() << " " << points_round_y->size() << points_round_z->size() << std::endl;
  /*return false;
  std::cout << (*points)[a1] << " " << (*points)[a2] << " " << (*points)[b1] << " " << (*points)[b2] << std::endl;
  std::cout << (*points_round_x)[a1] << " " << (*points_round_x)[a2] << " " << (*points_round_x)[b1] << " " << (*points_round_x)[b2] << std::endl;
  std::cout << (*points_round_y)[a1] << " " << (*points_round_y)[a2] << " " << (*points_round_y)[b1] << " " << (*points_round_y)[b2] << std::endl;
  std::cout << (*points_round_z)[a1] << " " << (*points_round_z)[a2] << " " << (*points_round_z)[b1] << " " << (*points_round_z)[b2] << std::endl;
  return false;*/
  std::array<Point_3, 4> p = {(*points)[a1],(*points)[a2],(*points)[b1],(*points)[b2]};
  std::array<Point_3, 4> px = {(*points_round_x)[a1],(*points_round_x)[a2],(*points_round_x)[b1],(*points_round_x)[b2]};
  std::array<Point_3, 4> py = {(*points_round_y)[a1],(*points_round_y)[a2],(*points_round_y)[b1],(*points_round_y)[b2]};
  std::array<Point_3, 4> pz = {(*points_round_z)[a1],(*points_round_z)[a2],(*points_round_z)[b1],(*points_round_z)[b2]};
//return false;
  //Using SOS, we can assume that the two segments are not coplanar at the beginning and after rounding in x (resp y, z) and before to making 
  //epsilon tend to zero.
  CGAL::Orientation orientation_begin=CGAL::orientation(p[0],p[1],p[2],p[3]);
  CGAL::Orientation orientation_round_x = CGAL::orientation(px[0],px[1],px[2],px[3]);
  CGAL::Orientation orientation_round_y = CGAL::orientation(py[0],py[1],py[2],py[3]);
  CGAL::Orientation orientation_round_z = CGAL::orientation(pz[0],pz[1],pz[2],pz[3]);

  //round_on_x
  if( orientation_begin != CGAL::COPLANAR &&
  	  orientation_round_x != CGAL::COPLANAR &&
  	  orientation_round_x != orientation_begin)
  {
	Segment_2 a_on_x(Point_2(p[0].y(),p[0].z()), Point_2(p[1].y(), p[1].z()));
  	Segment_2 b_on_x(Point_2(p[2].y(),p[2].z()), Point_2(p[3].y(), p[3].z()));
	if(seg_do_intersect(a_on_x,b_on_x,simplicialComplexOption))
	{
		intersectDuring=MotionStep::X_MOVING;
		return true;  
	}
  }

  //round on y
  if(orientation_round_y != CGAL::COPLANAR &&
     orientation_round_y != orientation_round_x &&
     (orientation_round_x != CGAL::COPLANAR || 
	 	(orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_round_y)))
  {
    Segment_2 ax_on_y(Point_2(px[0].x(),px[0].z()), Point_2(px[1].x(), px[1].z()));
    Segment_2 bx_on_y(Point_2(px[2].x(),px[2].z()), Point_2(px[3].x(), px[3].z()));
	if(seg_do_intersect(ax_on_y,bx_on_y,simplicialComplexOption))
	{
		intersectDuring=MotionStep::Y_MOVING;
		return true;
	}
  }

  //round on z
  if(orientation_round_z != CGAL::COPLANAR &&
  	 orientation_round_z != orientation_round_y &&
  	(orientation_round_y != CGAL::COPLANAR || 
	(orientation_round_x != CGAL::COPLANAR && orientation_round_x != orientation_round_z) || 
	(orientation_round_x == CGAL::COPLANAR && orientation_begin != CGAL::COPLANAR && orientation_begin != orientation_round_z)))
  {

	Segment_2 ay_on_z(Point_2(py[0].x(),py[0].y()), Point_2(py[1].x(), py[1].y()));
	Segment_2 by_on_z(Point_2(py[2].x(),py[2].y()), Point_2(py[3].x(), py[3].y()));
	if(seg_do_intersect(ay_on_z,by_on_z,simplicialComplexOption))
	{
		intersectDuring=MotionStep::Z_MOVING;
		return true;
	}
  }

  //round epsilon
  if(orientation_round_z == CGAL::COPLANAR &&
  	!( pz[2] == pz[0] || pz[2] == pz[1] || pz[3] == pz[1] || pz[3] == pz[0]))
  {
	if(seg_do_intersect(Segment_3(pz[0],pz[1]),Segment_3(pz[2],pz[3]),simplicialComplexOption))
	{
		Segment_2 ay_on_z(Point_2(py[0].x(),py[0].y()), Point_2(py[1].x(), py[1].y()));
		Segment_2 by_on_z(Point_2(py[2].x(),py[2].y()), Point_2(py[3].x(), py[3].y()));
		intersectDuring=MotionStep::COPLANAR;
		return true;
	}
  } //*/

  return false;
}

//_______________________________________________________________________________

/*
struct Functor_points_intersect{
	Functor_points_intersect(ListIntersection *list_intersections__):list_intersections_(list_intersections__){}
	
	void operator()(const Face_Vertex_Box& a_it, const Face_Vertex_Box& b_it) {
		//ListIntersection &list_intersections = *list_intersections_;
		num_tour++;
		
		P_index pa = a_it.handle()->vertex_parameter;
		P_index pb = b_it.handle()->vertex_parameter;
		if( points_do_intersect(points[pa], points[pb])){
			//list_intersections.face_vertex_intersections.emplace_back(b.second,pi,points[pi].slab(),points[pi].is_dummy, res_g);
			if(verbose){
				std::cout << "Intersection Vertex Vertex " << std::endl;
				std::cout << points[pa]() << " | " << points[pb]() << std::endl;
				std::cout << points[pa].round_to_integer() << " " << points[pb].round_to_integer() << std::endl << std::endl;
			}
			assert(!points[pa].round_to_integer() || !points[pb].round_to_integer());
			points[pa].set_round_to_integer();
			points[pb].set_round_to_integer();
		}
	}
	ListIntersection *list_intersections_;
	
};//*/


bool check_intersections(std::vector<Point_3> &points, std::vector< std::vector<size_t> > &triangles, FunctorSnap &snap, std::vector<Intersection> &out_intersection, bool simplicialComplexOption){
	std::cout << "Check Intersection" << std::endl;
	
	std::vector< Parameter_Face_Vertex > face_parameters;
	std::set< Parameter_Edge > edge_parameters;
	std::vector< Parameter_Face_Vertex > vertex_parameters;

	std::vector< std::array< Point_3, 4> > points_rounds;
	std::vector< Point_3 > points_round_x;
	std::vector< Point_3 > points_round_y;
	std::vector< Point_3 > points_round_z;
	for(size_t i=0; i<points.size(); ++i){
		auto &p = points[i];
		points_round_x.push_back(snap.snap_x(p, i));
		points_round_y.push_back(snap.snap_xy(p, i));
		points_round_z.push_back(snap.snap_xyz(p, i));
	}

	for(size_t i=0; i<points.size(); ++i){
		vertex_parameters.emplace_back(P_index(i));
	}
	for(size_t i=0; i<triangles.size(); ++i){
		face_parameters.emplace_back(F_index(i));
		assert(triangles[i].size()==3);
	}
	//TODO avoid doublon
	for(auto &t : triangles){
		if(t[0]<t[1])
			edge_parameters.emplace(P_index(t[0]), P_index(t[1]));
		else
			edge_parameters.emplace(P_index(t[1]), P_index(t[0]));
		if(t[0]<t[2])
			edge_parameters.emplace(P_index(t[0]), P_index(t[2]));
		else
			edge_parameters.emplace(P_index(t[2]), P_index(t[0]));
		if(t[1]<t[2])
			edge_parameters.emplace(P_index(t[1]), P_index(t[2]));
		else
			edge_parameters.emplace(P_index(t[2]), P_index(t[1]));
		//edge_parameters.emplace(P_index(t[1]), P_index(t[2]));
		//edge_parameters.emplace(P_index(t[2]), P_index(t[0]));
	}

	std::vector< Face_Vertex_Box > face_boxes;
	std::vector< Edge_Box > edge_boxes;
	std::vector< Face_Vertex_Box > vertex_boxes;
	
	for(auto it=face_parameters.begin(); it!=face_parameters.end(); ++it){
		face_boxes.emplace_back(bbox(triangles[it->is_face_parameter], points, snap), it);
	}
	for(auto it=vertex_parameters.begin(); it!=vertex_parameters.end(); ++it)
		vertex_boxes.emplace_back(extending_bbox(points[it->vertex_parameter], it->vertex_parameter, snap), it);
	
	for(auto it=edge_parameters.begin(); it!=edge_parameters.end(); ++it)
		edge_boxes.emplace_back(extending_bbox(points[it->first], it->first, snap)+extending_bbox(points[it->second], it->second, snap), it);
		
	Functor_edge_intersect functor1(&points, &triangles, &out_intersection, &points_round_x, &points_round_y, &points_round_z, simplicialComplexOption);
    CGAL::box_self_intersection_d( edge_boxes.begin(), edge_boxes.end(), functor1);
    Functor_face_intersect functor2(&points, &triangles, &out_intersection,  &points_round_x, &points_round_y, &points_round_z, simplicialComplexOption);
    CGAL::box_intersection_d( vertex_boxes.begin(), vertex_boxes.end(), face_boxes.begin(), face_boxes.end(), functor2);

    return out_intersection.size()!=0;
}


#endif