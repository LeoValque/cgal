// Copyright (c) 2024 ?? (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Leo Valque
//


#ifndef CGAL_POLYGON_MESH_PROCESSING_ROUNDING_VERTICES_H
#define CGAL_POLYGON_MESH_PROCESSING_ROUNDING_VERTICES_H

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <boost/core/bit.hpp>


#define CGAL_LINKED_WITH_TBB

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif

#define CGAL_PMP_ROUNDING_VERTICES_USE_DEFAULT_VERBOSE

#ifdef CGAL_PMP_ROUNDING_VERTICES_USE_DEFAULT_VERBOSE
#define CGAL_PMP_ROUNDING_VERTICES_VERBOSE(X) std::cout << X << "\n";
#endif

#ifndef CGAL_PMP_ROUNDING_VERTICES_VERBOSE
#define CGAL_PMP_ROUNDING_VERTICES_VERBOSE(MSG)
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

namespace internal{

    template <class K>
    typename K::Point_3 round_coordinates_to_double(const typename K::Point_3 &p){
        return typename K::Point_3(to_double(p.x()), to_double(p.y()), to_double(p.z()));
    }

    template <class K>
    double ceil(typename K::FT v){
        return std::ceil(CGAL::to_double(v));
    }

    template <class K>
    double exact_ceil(typename K::FT v){
        //Check if the precision of double is enough to describe correctly the integer
        if(std::ceil(to_interval(v).first) != std::ceil(to_interval(v).second)){
            //Compute exact
            //CGAL::exact(a);
            //Check again
            if(std::ceil(to_interval(v).first) != std::ceil(to_interval(v).second)){
                double int_lbrack =std::ceil(CGAL::to_interval(v).first);
                double int_rbrack =std::ceil(CGAL::to_interval(v).second);
                //Check all integer in the brack 
                for(double i=int_lbrack; i<=int_rbrack; i++){
                    if(v <= typename K::FT(i)){
                        return i;
                    }
                }
                assert(0);
            }
        }
        return std::ceil(CGAL::to_interval(v).first);
    }

    template <class K>
    double snap(typename K::FT v, double scale){
        return ceil<K>((v-0.5)*scale)/scale;
    }

    template <class K>
    typename K::Point_3 snap(const typename K::Point_3 &p, std::array<double, 3> scale){
        return typename K::Point_3(snap<K>(p.x(), scale[0]), snap<K>(p.y(), scale[1]), snap<K>(p.z(), scale[2]));
    }

    template <class K, class PointRange>
    void snap_nearby_vertices(PointRange &points, const std::set<size_t> indexes, const std::array<double, 3> scale){
	    if(indexes.size()==0)
		    return;

        //Round indexes vertices
        for(size_t i : indexes)
            points[i] = snap<K>(points[i],scale);

        //Copy and sort the indexes
        std::vector<typename K::Point_3> sorted_indexes;
        sorted_indexes.reserve(indexes.size());
        for(size_t i : indexes)
            sorted_indexes.emplace_back(points[i]);
        std::sort(sorted_indexes.begin(), sorted_indexes.end());

        //For all points, check if their snap versions are equal to one of the indexes
        #ifdef CGAL_LINKED_WITH_TBB
            if(true)
            //if (parallel_execution)
            {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, points.size()),
                            [&](const tbb::blocked_range<size_t>& r) {
                                for (size_t pi = r.begin(); pi != r.end(); ++pi){
                                    typename K::Point_3 snap_p=snap<K>(points[pi],scale);
                                    if(std::binary_search(sorted_indexes.begin(), sorted_indexes.end(), snap_p)){
                                        points[pi]=snap_p;
                                    }
                                }
                            }
                            );

            } else
        #endif
        for(typename K::Point_3 &p: points){
            typename K::Point_3 snap_p=snap<K>(p,scale);
                if(std::binary_search(sorted_indexes.begin(), sorted_indexes.end(), snap_p))
                    p=snap_p;
        }

    }

}

template <class K, class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
bool round_vertices_triangle_soup(PointRange& soup_points,
                              TriangleRange& soup_triangles,
                              const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    //size_t nb_iter = choose_parameter<Max_nb_rounding_iteration>(get_parameter(np, internal_np::max_nb_rounding_iteration));
    //size_t nb_iter = choose_parameter<Rounding_precision>(get_parameter(np, internal_np::rounding_precision));
    size_t nb_iter=20;
    size_t nb_bits=15;

    typedef std::vector< std::pair<std::size_t,std::size_t> > TriangleIdPairOutputVector;
    typedef std::size_t Input_TID;
    typedef std::size_t Input_VID;
    typedef std::pair<Input_TID, Input_TID> Pair_of_triangle_ids;

    //Compute the largest absolute value on each coordinate and take the smallest above power of two for each coordinate
    Bbox_3 bb=bbox_3(soup_points.begin(), soup_points.end());
    auto exponent = [](double v){
        int n;
        frexp(v, &n);
        return n;
    };
    std::array<double, 3> max = {std::max( std::abs(bb.xmin()),std::abs(bb.xmax())),
                                 std::max( std::abs(bb.ymin()),std::abs(bb.ymax())),
                                 std::max( std::abs(bb.zmin()),std::abs(bb.zmax()))};
    std::array<double, 3> scale = {std::pow(2.,nb_bits-exponent(max[0])),
                                   std::pow(2.,nb_bits-exponent(max[1])),
                                   std::pow(2.,nb_bits-exponent(max[2]))};

	for(size_t k=0; k<nb_iter; ++k)
    {
        //Round all coordinates on doubles
		for(typename K::Point_3 &p : soup_points)
			p=internal::round_coordinates_to_double<K>(p);

        //Compute pair of intersecting triangles
		std::vector<Pair_of_triangle_ids> si_pairs;
        repair_polygon_soup(soup_points, soup_triangles, np);
		triangle_soup_self_intersections(soup_points, soup_triangles, std::back_inserter(si_pairs), np);

        //If not intersection, it is terminate.
		if(si_pairs.empty())
        {
			return true;
		}

        // List all intersecting triangles and their vertices
		std::set<Input_VID> inter_points;
		std::set<Input_TID> inter_faces;
		for(Pair_of_triangle_ids pair: si_pairs)
        {
			inter_faces.emplace(pair.first);
			inter_faces.emplace(pair.second);
			for(size_t i=0; i<3; ++i){
				inter_points.emplace(soup_triangles[pair.first][i]);
				inter_points.emplace(soup_triangles[pair.second][i]);
			}
		}

        //Round coordinates of vertices nearby interecting ones to the integer grid
		internal::snap_nearby_vertices<K>(soup_points, inter_points, scale);

        //Refine self-intersections
		repair_polygon_soup(soup_points, soup_triangles, np);
		autorefine_triangle_soup(soup_points, soup_triangles, np);
	}

    //Fail to round without self-intersections
	return false;
}

}
}

#endif