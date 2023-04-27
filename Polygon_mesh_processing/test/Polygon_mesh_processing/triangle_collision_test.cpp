#define CGAL_NO_CDT_2_WARNING

#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <unistd.h>
//#include <CGAL/IO/STL_reader.h>
//#include "checker.hpp"
//#include "distance1.hpp"
#include <CGAL/Polygon_mesh_processing/triangle_collision.h>
#include <CGAL/Polygon_mesh_processing/proximity.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Aff_transformation_3.h>

#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <cmath>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
//#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
//#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Cartesian;
typedef Cartesian::Point_3 Double_Point_3;

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

//*
struct FunctorFloat : FunctorSnap{
    Kernel::FT snap_FT(Kernel::FT x, size_t pi){
        return round_to_float(x);
    }
};

struct FunctorDouble : FunctorSnap{
    Kernel::FT snap_FT(Kernel::FT x, size_t pi){
        return round_to_double(x);
    }
};//*/

void test_on_path(std::string path){
  namespace PMP = CGAL::Polygon_mesh_processing;

    std::vector< Double_Point_3 > double_points_triangle;
	std::vector< std::vector<size_t> > faces_triangles;
    std::cout << path  << std::endl;
  
	if(CGAL::IO::read_polygon_soup(path,double_points_triangle, faces_triangles)){
    PMP::repair_polygon_soup(double_points_triangle, faces_triangles);
    CGAL::Surface_mesh<Point_3> mesh;
    PMP::orient_polygon_soup(double_points_triangle, faces_triangles);
    //CGAL::Cartesian_converter_property_map<Double_Point_3, decltype(mesh.points())> vpm(mesh.points());
    auto vpm=CGAL::make_cartesian_converter_property_map<Double_Point_3>(mesh.points());
    PMP::polygon_soup_to_polygon_mesh(double_points_triangle, faces_triangles, mesh, CGAL::parameters::all_default(), CGAL::parameters::vertex_point_map(vpm));
    PMP::triangulate_faces(mesh);
    if(PMP::does_self_intersect(mesh)){
      try{
        //PMP::autorefine(mesh);
      } catch(PMP::Corefinement::Triple_intersection_exception&()){
        //DEBUG_OUT;
        exit(1);
        //PMP::autorefine(mesh);
      }
    }

    faces_triangles.clear();
    std::vector< Point_3 > points_triangle;
    PMP::polygon_mesh_to_polygon_soup(mesh, points_triangle, faces_triangles);
    //exit(1);
  //*
    std::vector<Intersection> out_intersection;
    FunctorFloat snap;
    FunctorDouble snap2;
	  check_intersections(points_triangle, faces_triangles, snap, out_intersection, true);
    //check_intersections(points_triangle, faces_triangles, snap, out_intersection, true);


    std::cout << "Number of Intersection: " << out_intersection.size() << std::endl;
    for(auto &intersection : out_intersection){
      if(intersection.is_edge_edge){
        std::cout << "Intersection between edges: " << intersection.ei1.first << "->" << intersection.ei1.second << " , " << intersection.ei2.first << "->" << intersection.ei2.second << " during " << static_cast<int> (intersection.intersectionDuring) << std::endl;
        std::cout << points_triangle[intersection.ei1.first] << " -> " << points_triangle[intersection.ei1.second] << std::endl;
        std::cout << points_triangle[intersection.ei2.first] << " -> " << points_triangle[intersection.ei2.second] << "\n" << std::endl;
      } else {
        std::cout << "Intersection between a face and a vertex: " << intersection.pi << " , " << intersection.fi << ": "<< faces_triangles[intersection.fi][0] << "->" << faces_triangles[intersection.fi][1] << "->" << faces_triangles[intersection.fi][2] << " during " << static_cast<int>(intersection.intersectionDuring) << std::endl;
        std::cout << points_triangle[faces_triangles[intersection.fi][0]] << " -> " << points_triangle[faces_triangles[intersection.fi][1]] << " -> " << points_triangle[faces_triangles[intersection.fi][2]] << std::endl;
        std::cout << points_triangle[intersection.pi] << "\n" << std::endl;
      }
    }//*/

	} else{
		std::cerr << "Not found: " << path << std::endl;
	}
}

int main(int argc, char *argv[]){
  std::cout.precision(17);

  //std::cout << argc << " " << argv[1] << " " << (argv[1]=="test_thingi") << " " << argv[2] << std::endl;
  if(argc==3){
  	//DEBUG_OUT;
    if(std::string(argv[1])==std::string("run_on_thingi"))
  	  test_on_path("../../data/Thingi10K/raw_meshes/"+std::string(argv[2]));
    else
      test_on_path(std::string(argv[2]));
  	return 0;
  }
  return 0;
}
