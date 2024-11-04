
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/rounding_vertices.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Real_timer.h>

#include <boost/container/small_vector.hpp>

#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  std::vector<Point> input_points;
  std::vector<boost::container::small_vector<std::size_t, 3>> input_triangles;
  if (!CGAL::IO::read_polygon_soup(filename, input_points, input_triangles))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }
  PMP::repair_polygon_soup(input_points, input_triangles);
  PMP::triangulate_polygons(input_points, input_triangles);

  CGAL::Real_timer t;
  t.start();
  PMP::round_vertices_triangle_soup(input_points, input_triangles,
                                //CGAL::parameters::concurrency_tag(CGAL::Sequential_tag()));
                                //CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
                                CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));
  t.stop();
  std::cout << "#points = " << input_points.size() << " and #triangles = " << input_triangles.size() << " in " << t.time() << " sec." << std::endl;

  std::vector<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> output_points;
  for(auto &p: input_points)
    output_points.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
  CGAL::IO::write_polygon_soup("autorefined_and_rounded.off", output_points, input_triangles, CGAL::parameters::stream_precision(17));

  return 0;
}
