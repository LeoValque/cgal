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

#ifndef CGAL_POLYGON_MESH_PROCESSING_PROXIMITY_H
#define CGAL_POLYGON_MESH_PROCESSING_PROXIMITY_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <filesystem>
#include <limits>
#include <CGAL/box_intersection_d.h>

#include <filesystem>
#include <limits>

#define DEBUG_OUT std::cout << "Debug: " << std::filesystem::path(__FILE__).filename().string() << ":" << __LINE__ << " " << __func__ <<std::endl

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

bool at_distance_one(const Triangle_3 a, const Triangle_3 b);
bool triangle_has_on_projected_point(const Triangle_3 &tr, const Point_3 &p);

double round_to_float(Kernel::FT a);
double round_to_double(Kernel::FT a);

void distance1features(std::vector<Point_3> &points, std::vector< std::vector<size_t> > &triangles, double r, std::vector< std::pair<size_t, size_t> > &outs, bool considerAdjacentFeature);

typedef std::vector< size_t >::iterator Iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;

Bbox_3 extend(Bbox_3 b, double a){
	return Bbox_3(b.xmin(), b.ymin(), b.zmin(), b.xmax()+a, b.ymax()+a, b.zmax()+a);
}

bool at_distance_one(const Triangle_3 a, const Triangle_3 b, const double r){
    //Kernel::FT three(3);
    //Kernel::Compare_squared_distance_3 compare_squared_distance;
    //return CGAL::possibly(compare_squared_distance(a,b,three)==CGAL::SMALLER);


    Kernel::Compare_squared_distance_3 compare_squared_distance;
    Kernel::FT three_r(3*r);
	//First, we check the distance between each pair of segments
	//After, we check the distance between each vertew to the opposite triangle
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j)
            //If the distance between the lines is more than 3, the distance between the segments is more than 3 and thus, we not compute it.
            //If the distance between the lines is less than 3, it is computed again in the distance between the segments but this gains a bit of time because 
            //most pair of lines are at distance more than 3
            if(CGAL::possibly(compare_squared_distance(Line_3(a[i],a[(i+1)%3]), Line_3(b[j],b[(j+1)%3]),three_r)==CGAL::SMALLER) && 
               CGAL::possibly(compare_squared_distance(Segment_3(a[i],a[(i+1)%3]), Segment_3(b[j],b[(j+1)%3]),three_r)==CGAL::SMALLER))
            	return true;

        //If the distance to the plane is more than 3, the distance to the triangle is more than 3 and this filters many cases.
        //We already test that all edges are at distance more than 3 so we can simply test if the projection of the vertex is inside the triangle.
        if((CGAL::possibly(compare_squared_distance(a.supporting_plane(),b[i],three_r)==CGAL::SMALLER) && triangle_has_on_projected_point(a,b[i])==CGAL::SMALLER)
           || (CGAL::possibly(compare_squared_distance(a[i],b.supporting_plane(),three_r)==CGAL::SMALLER) && triangle_has_on_projected_point(b,a[i])==CGAL::SMALLER))
    		return true;
    }
    return false;
}

class FunctorDistance1{
public:
    FunctorDistance1(std::vector<Point_3> *points_, std::vector< std::vector<size_t> > *triangles_, double r_, std::vector< std::pair<size_t, size_t> > *outs_, bool cAdja_):
		points(points_), triangles(triangles_), r(r_), outs(outs_), cAdja(cAdja_){}
    void operator()(const Box& a, const Box& b) {
		std::vector<size_t> &vector_a = (*triangles)[*(a.handle())];
		std::vector<size_t> &vector_b = (*triangles)[*(b.handle())];
		for(auto i : vector_a)
			for(auto j : vector_b)
				if(i==j){
					if(cAdja)
						outs->emplace_back(*(a.handle()), *(b.handle()));
					return;
				}
		if(at_distance_one(Triangle_3((*points)[vector_a[0]], (*points)[vector_a[1]], (*points)[vector_a[2]]), Triangle_3((*points)[vector_b[0]], (*points)[vector_b[1]], (*points)[vector_b[2]]), r))
			outs->emplace_back(*(a.handle()), *(b.handle()));
    }

	std::vector<Point_3> *points;
	std::vector< std::vector<size_t> > *triangles;
	double r; 
	std::vector< std::pair<size_t, size_t> > *outs;
	bool cAdja;
};

void proximity_triangle_soup(std::vector<Point_3> &points, std::vector< std::vector<size_t> > &triangles, double r, std::vector< std::pair<size_t, size_t> > &outs, bool considerAdjacentFeature){
    std::vector<size_t> triangleIDs(triangles.size());
	for(size_t i=0; i!=triangles.size(); ++i)
		triangleIDs[i]=i;
	std::vector< Box > faceBoxes;
    for(Iterator it = triangleIDs.begin(); it != triangleIDs.end(); ++it)
        faceBoxes.emplace_back(extend(points[triangles[(*it)][0]].bbox()+points[triangles[(*it)][1]].bbox()+points[triangles[(*it)][2]].bbox(),r), it);

    FunctorDistance1 functor(&points, &triangles, r, &outs, considerAdjacentFeature);
    CGAL::box_self_intersection_d( faceBoxes.begin(), faceBoxes.end(), functor);
}

bool triangle_has_on_projected_point(const Triangle_3 &tr, const Point_3 &p){
    Point_3 p_project = tr.supporting_plane().projection(p);
    return CGAL::possibly(tr.has_on(p_project));
}

//This function return true if two triangles are at L1-distance less than 1. Being at L2-distance less than sqrt(3) implies the former.
//possibly avoid an exact computation and may return true two triangles are at distance close and above to 1.
bool at_distance_one(const Triangle_3 a, const Triangle_3 b){
    //Kernel::FT three(3);
    //Kernel::Compare_squared_distance_3 compare_squared_distance;
    //return CGAL::possibly(compare_squared_distance(a,b,three)==CGAL::SMALLER);


    Kernel::Compare_squared_distance_3 compare_squared_distance;
    Kernel::FT three(3);
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j)
            //If the distance between the lines is more than 3, the distance between the segments is more than 3 and thus, we not compute it.
            //If the distance between the lines is less than 3, it is computed again in the distance between the segments but this gains a bit of time because 
            //most pair of lines are at distance less than 3
            if(CGAL::possibly(compare_squared_distance(Line_3(a[i],a[(i+1)%3]), Line_3(b[j],b[(j+1)%3]),three)==CGAL::SMALLER) && 
               CGAL::possibly(compare_squared_distance(Segment_3(a[i],a[(i+1)%3]), Segment_3(b[j],b[(j+1)%3]),three)==CGAL::SMALLER))
            	return true;

        //If the distance to the plane is more than 3, the distance to the triangle is more than 3 and this filters many cases.
        //We already test that all edges are at distance more than 3 so we can simply test if the vertex is inside the triangle.
        if((CGAL::possibly(compare_squared_distance(a.supporting_plane(),b[i],three)==CGAL::SMALLER) && triangle_has_on_projected_point(a,b[i]))
           || (CGAL::possibly(compare_squared_distance(a[i],b.supporting_plane(),three)==CGAL::SMALLER) && triangle_has_on_projected_point(b,a[i])))
    		return true;
    }
    return false;
}

double round_to_double(Kernel::FT a){
    a.set_relative_precision_of_to_double(1/std::pow(2,52));
    //a.set_relative_precision_of_to_double(1/std::pow(2,26));
    return CGAL::to_double(a);
}

double round_to_float(Kernel::FT a){
    a.set_relative_precision_of_to_double(1/std::pow(2,26));
    //a.set_relative_precision_of_to_double(1/std::pow(2,26));
    return (float) CGAL::to_double(a);
}

#endif