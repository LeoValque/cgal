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
#include <CGAL/box_intersection_d.h>

typedef CGAL::Bbox_2       Bbox_2;
typedef CGAL::Bbox_3       Bbox_3;

typedef std::vector< size_t >::iterator Iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;

Bbox_3 extend(Bbox_3 b, double a){
	return Bbox_3(b.xmin(), b.ymin(), b.zmin(), b.xmax()+a, b.ymax()+a, b.zmax()+a);
}

template<class Kernel>
bool triangle_has_on_projected_point(const typename Kernel::Triangle_3 &tr, const typename Kernel::Point_3 &p){
    typename Kernel::Point_3 p_project = tr.supporting_plane().projection(p);
    return CGAL::possibly(tr.has_on(p_project));
}

template <class Kernel>
bool at_distance_one_bis(const typename Kernel::Triangle_3 a, const typename Kernel::Triangle_3 b, const double r){
    typename Kernel::FT three_r(3*r);
    typename Kernel::Compare_squared_distance_3 compare_squared_distance;
    return CGAL::possibly(compare_squared_distance(a,b,three_r)==CGAL::SMALLER);
}

//This function return true if two triangles are at L1-distance less than r. Being at L2-distance less than sqrt(3)*r implies the former.
//possibly avoid an exact computation and may return true two triangles are at distance close and above to r.
template <class Kernel>
bool at_distance_one(const typename Kernel::Triangle_3 a, const typename Kernel::Triangle_3 b, const double r){
    typename Kernel::Compare_squared_distance_3 compare_squared_distance;
    typename Kernel::FT three_r(3*r);

	//First, we check the distance between each pair of segments
	//After, we check the distance between each vertew to the opposite triangle
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j)
            //If the distance between the lines is more than 3, the distance between the segments is more than 3 and thus, we not compute it.
            //If the distance between the lines is less than 3, it is computed again in the distance between the segments but this gains a bit of time because 
            //most pair of lines are at distance more than 3
            if(CGAL::possibly(compare_squared_distance(typename Kernel::Line_3(a[i],a[(i+1)%3]), typename Kernel::Line_3(b[j],b[(j+1)%3]),three_r)==CGAL::SMALLER) && 
               CGAL::possibly(compare_squared_distance(typename Kernel::Segment_3(a[i],a[(i+1)%3]), typename Kernel::Segment_3(b[j],b[(j+1)%3]),three_r)==CGAL::SMALLER))
            	return true;

        //If the distance to the plane is more than 3, the distance to the triangle is more than 3 and this filters many cases.
        //We already test that all edges are at distance more than 3 so we can simply test if the projection of the vertex is inside the triangle.
        if((CGAL::possibly(compare_squared_distance(a.supporting_plane(),b[i],three_r)==CGAL::SMALLER) && triangle_has_on_projected_point<Kernel>(a,b[i]))
           || (CGAL::possibly(compare_squared_distance(a[i],b.supporting_plane(),three_r)==CGAL::SMALLER) && triangle_has_on_projected_point<Kernel>(b,a[i])))
    		return true;
    }
    return false;
}

template< class Kernel>
class FunctorDistance1Bis{
public:
    FunctorDistance1Bis(const std::vector< typename Kernel::Point_3> *points_, const std::vector< std::vector<size_t> > *triangles_, const double r_, std::vector< std::pair<size_t, size_t> > *outs_, bool cAdja_):
		points(points_), triangles(triangles_), r(r_), outs(outs_), cAdja(cAdja_){}
    void operator()(const Box& a, const Box& b) {
		const std::vector<size_t> &vector_a = (*triangles)[*(a.handle())];
		const std::vector<size_t> &vector_b = (*triangles)[*(b.handle())];
		for(auto i : vector_a)
			for(auto j : vector_b)
				if(i==j){
					if(cAdja)
						outs->emplace_back(*(a.handle()), *(b.handle()));
					return;
				}
		if(at_distance_one_bis<Kernel>(typename Kernel::Triangle_3((*points)[vector_a[0]], (*points)[vector_a[1]], (*points)[vector_a[2]]), 
                               typename Kernel::Triangle_3((*points)[vector_b[0]], (*points)[vector_b[1]], (*points)[vector_b[2]]), r))
			outs->emplace_back(*(a.handle()), *(b.handle()));
    }

	const std::vector<typename Kernel::Point_3> *points;
	const std::vector< std::vector<size_t> > *triangles;
	double r; 
	std::vector< std::pair<size_t, size_t> > *outs;
	bool cAdja;
};

template<class Kernel>
class FunctorDistance1{
public:
    FunctorDistance1(const std::vector<typename Kernel::Point_3> *points_, const std::vector< std::vector<size_t> > *triangles_, double r_, std::vector< std::pair<size_t, size_t> > *outs_, bool cAdja_):
		points(points_), triangles(triangles_), r(r_), outs(outs_), cAdja(cAdja_){}
    void operator()(const Box& a, const Box& b) {
		const std::vector<size_t> &vector_a = (*triangles)[*(a.handle())];
		const std::vector<size_t> &vector_b = (*triangles)[*(b.handle())];
		for(auto i : vector_a)
			for(auto j : vector_b)
				if(i==j){
					if(cAdja)
						outs->emplace_back(*(a.handle()), *(b.handle()));
					return;
				}
		if(at_distance_one<Kernel>(typename Kernel::Triangle_3((*points)[vector_a[0]], (*points)[vector_a[1]], (*points)[vector_a[2]]), 
                           typename Kernel::Triangle_3((*points)[vector_b[0]], (*points)[vector_b[1]], (*points)[vector_b[2]]), r))
			outs->emplace_back(*(a.handle()), *(b.handle()));
    }

	const std::vector<typename Kernel::Point_3> *points;
	const std::vector< std::vector<size_t> > *triangles;
	double r; 
	std::vector< std::pair<size_t, size_t> > *outs;
	bool cAdja;
};

template<class Kernel>
void proximity_triangle_soup(const std::vector<typename Kernel::Point_3> &points, const std::vector< std::vector<size_t> > &triangles, const double r, std::vector< std::pair<size_t, size_t> > &outs, bool considerAdjacentFeature){
    std::vector<size_t> triangleIDs(triangles.size());
	for(size_t i=0; i!=triangles.size(); ++i)
		triangleIDs[i]=i;
	std::vector< Box > faceBoxes;
    for(Iterator it = triangleIDs.begin(); it != triangleIDs.end(); ++it)
        faceBoxes.emplace_back(extend(points[triangles[(*it)][0]].bbox()+points[triangles[(*it)][1]].bbox()+points[triangles[(*it)][2]].bbox(),r), it);

    FunctorDistance1<Kernel> functor(&points, &triangles, r, &outs, considerAdjacentFeature);
    CGAL::box_self_intersection_d( faceBoxes.begin(), faceBoxes.end(), functor);
}

template<class Kernel>
void proximity_triangle_soup_bis(std::vector<typename Kernel::Point_3> &points, std::vector< std::vector<size_t> > &triangles, double r, std::vector< std::pair<size_t, size_t> > &outs, bool considerAdjacentFeature){
    std::vector<size_t> triangleIDs(triangles.size());
	for(size_t i=0; i!=triangles.size(); ++i)
		triangleIDs[i]=i;
	std::vector< Box > faceBoxes;
    for(Iterator it = triangleIDs.begin(); it != triangleIDs.end(); ++it)
        faceBoxes.emplace_back(extend(points[triangles[(*it)][0]].bbox()+points[triangles[(*it)][1]].bbox()+points[triangles[(*it)][2]].bbox(),r), it);

    FunctorDistance1Bis<Kernel> functor(&points, &triangles, r, &outs, considerAdjacentFeature);
    CGAL::box_self_intersection_d( faceBoxes.begin(), faceBoxes.end(), functor);
}

#endif