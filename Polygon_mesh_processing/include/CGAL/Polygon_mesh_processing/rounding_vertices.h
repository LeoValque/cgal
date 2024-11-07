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
#include <CGAL/Lazy_exact_nt.h>
#include <boost/core/bit.hpp>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#endif

#ifdef CGAL_PMP_ROUNDING_VERTICES_USE_DEFAULT_VERBOSE
#define CGAL_PMP_ROUNDING_VERTICES_VERBOSE(X) std::cout << X << "\n";
#endif

#ifndef CGAL_PMP_ROUNDING_VERTICES_VERBOSE
#define CGAL_PMP_ROUNDING_VERTICES_VERBOSE(MSG)
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

namespace internal{

    template <class NT>
    double ceil(NT v){
        return std::ceil(to_double(v));
    }

    //If Lazy_exact is not use, they are no difference between ceil/exact_ceil and to_double/to_exact_closest_double
    template <class NT>
    double to_exact_closest_double(NT v){
        return to_double(v);
    }

    template <class NT>
    double exact_ceil(NT v){
        return ceil(v);
    }

    //If Lazy_exact is use, we can filter if the interval is large of one ulp
    template <class NT>
    double to_exact_closest_double(Lazy_exact_nt< NT > v){
        if( std::nextafter(to_interval(v).first, to_interval(v).second)==to_interval(v).second )
            return to_interval(v).first;
        return to_double(exact(v));
    }
    
    //If Lazy_exact is use, we can filter with the interval if they have same ceil values
    template <class NT>
    double exact_ceil(Lazy_exact_nt< NT > v){
        //Check if the precision of double is enough to describe correctly the integer
        if(std::ceil(to_interval(v).first) != std::ceil(to_interval(v).second)){
            //Compute exact and check again
            exact(v);
            if(std::ceil(to_interval(v).first) != std::ceil(to_interval(v).second)){
                //Test all integers in the bracket in increasing order
                double int_lbrack =std::ceil(CGAL::to_interval(v).first);
                double int_rbrack =std::ceil(CGAL::to_interval(v).second);
                for(double i=int_lbrack; i<=int_rbrack; i++){
                    if(v <= i){
                        return i;
                    }
                }
                assert(0);
            }
        }
        return std::ceil(CGAL::to_interval(v).first);
    }

}

/**
*
* Check if the coordinates of all the vertices fit in double. 
*
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param soup_points points of the soup of polygons
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
*     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
*     \cgalParamDefault{`CGAL::Sequential_tag`}
*   \cgalParamNEnd
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `soup_points`}
*     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return true if all the coordinates fit in double, false if not
*/

template <class PointRange, class NamedParameters = parameters::Default_named_parameters>
bool does_triangle_soup_fit_in_double(PointRange& soup_points,
                                 const NamedParameters& np = parameters::default_values())
{
    typedef typename GetPolygonSoupGeomTraits<PointRange, NamedParameters>::type K;

    double v;
    for(typename K::Point_3 &p: soup_points)
         if( (to_double(p.x())!=p.x()) || (to_double(p.y())!=p.y()) && (to_double(p.z())!=p.z()))
            return false;
    return true;
}

/**
*
* Refines a soup of triangles and round the coordinates of the vertices so that no pair of triangles intersects and all the coordinates fit in double. 
* The function iterates as long as self-intersection occurs or the max number of iteration is reached.
* Output triangles may share a common edge or a common vertex (but with the same indexed position in `points`).
* Note that the function calls repair_polygon_soup and thus the point and polygon containers will be modified by the repairing operations, and thus the indexing of the polygons will also be changed.
* Note that the output is not guarantee to be intersection free.
*
* @tparam PointRange a model of the concept `RandomAccessContainer`
* whose value type is the point type
* @tparam TriangleRange a model of the concepts `RandomAccessContainer`, `BackInsertionSequence` and `Swappable`, whose
* value type is a model of the concept `RandomAccessContainer` whose value type is convertible to `std::size_t` and that
* is constructible from an `std::initializer_list<std::size_t>` of size 3.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param soup_points points of the soup of polygons
* @param soup_triangles each element in the range describes a triangle using the indexed position of the points in `soup_points`
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{concurrency_tag}
*     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
*     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
*     \cgalParamDefault{`CGAL::Sequential_tag`}
*   \cgalParamNEnd
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `soup_points`}
*     \cgalParamType{a model of `ReadWritePropertyMap` whose value type is a point type}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the point type.}
*   \cgalParamNEnd
*   \cgalParamNBegin{max_nb_iteration}
*     \cgalParamDescription{The maximum number of iteration of the algorithm.}
*     \cgalParamType{size_t}
*     \cgalParamDefault{20}
*   \cgalParamNEnd
*   \cgalParamNBegin{rounding_precision}
*     \cgalParamDescription{Number of digits of the coordinates conserved during each iteration of the rounding process}
*     \cgalParamType{size_t}
*     \cgalParamDefault{23}
*     \cgalParamExtra{The value must be lower than 52 (number of digits of a mantissa). A higher value reduce the distance between a point and its rounding value but increase the chance of the output to be self-intersecting}
*   \cgalParamNEnd
* \cgalParamNBegin{exact_rounding}
*     \cgalParamDescription{Compute the exact value before to round a coordinate giving guarantee on the distance between a point and its rounding value.}
*     \cgalParamType{bool}
*     \cgalParamDefault{true}
*     \cgalParamExtra{Use only if the number type is a Lazy_exact_NT}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @return true if the output is intersection free, false if not
*/

template <class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
bool round_vertices_triangle_soup(PointRange& soup_points,
                              TriangleRange& soup_triangles,
                              const NamedParameters& np = parameters::default_values())
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef std::vector< std::pair<std::size_t,std::size_t> > TriangleIdPairOutputVector;
    typedef std::size_t Input_TID;
    typedef std::size_t Input_VID;
    typedef std::pair<Input_TID, Input_TID> Pair_of_triangle_ids;

    typedef typename GetPolygonSoupGeomTraits<PointRange, NamedParameters>::type K;
    typedef typename internal_np::Lookup_named_param_def <
        internal_np::concurrency_tag_t,
        NamedParameters,
        Sequential_tag
    > ::type Concurrency_tag;
    
    constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

    //TODO
    //constexpr size_t nb_iter = choose_parameter<Max_nb_rounding_iteration>(get_parameter(np, internal_np::max_nb_rounding_iteration));
    //constexpr size_t nb_bits = choose_parameter<Rounding_precision>(get_parameter(np, internal_np::rounding_precision));
    //constexpr bool exact_rounding = choose_parameter<Exact_rounding>(get_parameter(np, internal_np::rounding_precision));
    constexpr size_t nb_iter=20;
    constexpr size_t nb_bits=23;
    constexpr bool exact_rounding=true;

    //Compute the largest absolute value on each coordinate and take the smallest above power of two for each coordinate
    CGAL_PMP_ROUNDING_VERTICES_VERBOSE("compute scaling of the coordinates")
    Bbox_3 bb=bbox_3(soup_points.begin(), soup_points.end());
    auto exponent = [](const double v){
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

    //Define snap functions
    auto snap_v=[exact_rounding](typename K::FT v, double scale){
        if(exact_rounding)
            return internal::exact_ceil((v-0.5)*scale)/scale;
        else
            return internal::ceil((v-0.5)*scale)/scale;
    };

    auto snap=[&scale, &snap_v](const typename K::Point_3 &p){
        return typename K::Point_3(snap_v(p.x(), scale[0]), snap_v(p.y(), scale[1]), snap_v(p.z(), scale[2]));
    };

    CGAL_PMP_ROUNDING_VERTICES_VERBOSE("start the while loop")
	for(size_t k=0; k<nb_iter; ++k)
    {
        CGAL_PMP_ROUNDING_VERTICES_VERBOSE("start iteration " << (k+1))

        //Round all coordinates on doubles
        CGAL_PMP_ROUNDING_VERTICES_VERBOSE("Round coordinates")
        if(exact_rounding)
        {
            for(typename K::Point_3 &p : soup_points)
		        p=typename K::Point_3(internal::to_exact_closest_double(p.x()), 
                                      internal::to_exact_closest_double(p.y()), 
                                      internal::to_exact_closest_double(p.z()));
        } else {
		    for(typename K::Point_3 &p : soup_points)
			    p=typename K::Point_3(to_double(p.x()), to_double(p.y()), to_double(p.z()));
        }

        //Compute pair of intersecting triangles
        CGAL_PMP_ROUNDING_VERTICES_VERBOSE("Compute pairs of intersecting triangles")
		std::vector<Pair_of_triangle_ids> si_pairs;
        repair_polygon_soup(soup_points, soup_triangles, np);
		triangle_soup_self_intersections(soup_points, soup_triangles, std::back_inserter(si_pairs), np);

        //If not intersection, it is terminate.
		if(si_pairs.empty())
        {
            #ifndef CGAL_NDEBUG
            CGAL_PMP_ROUNDING_VERTICES_VERBOSE("check soup");
            CGAL_assertion( does_triangle_soup_fit_in_double<K>(soup_points, soup_triangles) );
            CGAL_assertion( !does_triangle_soup_self_intersect<Concurrency_tag>(soup_points, soup_triangles) );
            #endif
            CGAL_PMP_ROUNDING_VERTICES_VERBOSE("Done")
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

        //Snap the vertices of intersecting triangles
        CGAL_PMP_ROUNDING_VERTICES_VERBOSE("Snap the vertices of intersecting triangles")
        for(size_t i : inter_points)
            soup_points[i] = snap(soup_points[i]);

        //Copy and sort the snap vertices
        std::vector<typename K::Point_3> sorted_snap_vertices;
        sorted_snap_vertices.reserve(inter_points.size());
        for(size_t i : inter_points)
            sorted_snap_vertices.emplace_back(soup_points[i]);
        std::sort(sorted_snap_vertices.begin(), sorted_snap_vertices.end());

        //Snap the vertices nearby than these already snapped
        //For all the points, we check if their snap version are equal to one of the indexes
        CGAL_PMP_ROUNDING_VERTICES_VERBOSE("Snap vertices nearby than these already rounded")
        #ifdef CGAL_LINKED_WITH_TBB
            if(true)
            //if (parallel_execution)
            {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, soup_points.size()),
                            [&](const tbb::blocked_range<size_t>& r) {
                                for (size_t pi = r.begin(); pi != r.end(); ++pi){
                                    typename K::Point_3 snap_p=snap(soup_points[pi]);
                                    if(std::binary_search(sorted_snap_vertices.begin(), sorted_snap_vertices.end(), snap_p)){
                                        soup_points[pi]=snap_p;
                                    }
                                }
                            }
                            );

            } else
        #endif
        for(typename K::Point_3 &p: soup_points){
            typename K::Point_3 snap_p=snap(p);
                if(std::binary_search(sorted_snap_vertices.begin(), sorted_snap_vertices.end(), snap_p))
                    p=snap_p;
        }

        //Refine self-intersections
        CGAL_PMP_ROUNDING_VERTICES_VERBOSE("Refine self-intersections")
		repair_polygon_soup(soup_points, soup_triangles, np);
		autorefine_triangle_soup(soup_points, soup_triangles, np);
	}

    //Fail to round without self-intersections
	return false;
}

}
}

#endif