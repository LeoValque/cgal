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

namespace CGAL {
namespace Polygon_mesh_processing{

typedef std::vector< size_t >::iterator VectorIterator;
typedef Box_intersection_d::Box_with_handle_d<double,3,VectorIterator> Box;

template<class Kernel>
class FunctorDistance1{
public:
    FunctorDistance1(const std::vector<typename Kernel::Point_3> *points_, const std::vector< std::vector<size_t> > *triangles_, const typename Kernel::FT r_, std::vector< std::pair<size_t, size_t> > *outs_, bool cAdja_):
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

		typename Kernel::Compare_squared_distance_3 compare_squared_distance;
		if(compare_squared_distance(typename Kernel::Triangle_3((*points)[vector_a[0]], (*points)[vector_a[1]], (*points)[vector_a[2]]), 
                           typename Kernel::Triangle_3((*points)[vector_b[0]], (*points)[vector_b[1]], (*points)[vector_b[2]]), r)==CGAL::SMALLER)
			outs->emplace_back(*(a.handle()), *(b.handle()));
    }

	const std::vector<typename Kernel::Point_3> *points;
	const std::vector< std::vector<size_t> > *triangles;
	typename Kernel::FT r; 
	std::vector< std::pair<size_t, size_t> > *outs;
	bool cAdja;
};

template<class Kernel>
void proximity_pairs_in_triangles_soup(const std::vector<typename Kernel::Point_3> &points, const std::vector< std::vector<size_t> > &triangles, const typename Kernel::FT r, std::vector< std::pair<size_t, size_t> > &outs, bool considerAdjacentFeature, bool soup_self_intersect){
    auto extend=[](Bbox_3 b, double a){
		return Bbox_3(b.xmin(), b.ymin(), b.zmin(), b.xmax()+a, b.ymax()+a, b.zmax()+a);
	};
	
	std::vector<size_t> triangleIDs(triangles.size());
	for(size_t i=0; i!=triangles.size(); ++i)
		triangleIDs[i]=i;
	std::vector< Box > faceBoxes;
    for(VectorIterator it = triangleIDs.begin(); it != triangleIDs.end(); ++it)
        faceBoxes.emplace_back(extend(points[triangles[(*it)][0]].bbox()+points[triangles[(*it)][1]].bbox()+points[triangles[(*it)][2]].bbox(),CGAL::to_interval(r).second), it);

    FunctorDistance1<Kernel> functor(&points, &triangles, r, &outs, considerAdjacentFeature);
    box_self_intersection_d( faceBoxes.begin(), faceBoxes.end(), functor);
}





#if 0



template< class Kernel>
class FunctorDistance1Old{
public:
    FunctorDistance1Old(const std::vector< typename Kernel::Point_3> *points_, const std::vector< std::vector<size_t> > *triangles_, const typename Kernel::FT r_, std::vector< std::pair<size_t, size_t> > *outs_, bool cAdja_):
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
		if(CGAL::possibly(CGAL::compare( CGAL::squared_distance(typename Kernel::Triangle_3((*points)[vector_a[0]], (*points)[vector_a[1]], (*points)[vector_a[2]]), 
                               typename Kernel::Triangle_3((*points)[vector_b[0]], (*points)[vector_b[1]], (*points)[vector_b[2]])), r)==CGAL::SMALLER))
			outs->emplace_back(*(a.handle()), *(b.handle()));
    }

	const std::vector<typename Kernel::Point_3> *points;
	const std::vector< std::vector<size_t> > *triangles;
	typename Kernel::FT r; 
	std::vector< std::pair<size_t, size_t> > *outs;
	bool cAdja;
};

template<class Kernel>
void proximity_triangle_soup_old(const std::vector<typename Kernel::Point_3> &points, const std::vector< std::vector<size_t> > &triangles, const typename Kernel::FT r, std::vector< std::pair<size_t, size_t> > &outs, bool considerAdjacentFeature){
    auto extend=[](Bbox_3 b, double a){
		return Bbox_3(b.xmin(), b.ymin(), b.zmin(), b.xmax()+a, b.ymax()+a, b.zmax()+a);
	};
	
	std::vector<size_t> triangleIDs(triangles.size());
	for(size_t i=0; i!=triangles.size(); ++i)
		triangleIDs[i]=i;
	std::vector< Box > faceBoxes;
    for(VectorIterator it = triangleIDs.begin(); it != triangleIDs.end(); ++it)
        faceBoxes.emplace_back(extend(points[triangles[(*it)][0]].bbox()+points[triangles[(*it)][1]].bbox()+points[triangles[(*it)][2]].bbox(),CGAL::to_interval(r).second), it);

    FunctorDistance1Old<Kernel> functor(&points, &triangles, r, &outs, considerAdjacentFeature);
    CGAL::box_self_intersection_d( faceBoxes.begin(), faceBoxes.end(), functor);
}

#endif

}
}

#endif