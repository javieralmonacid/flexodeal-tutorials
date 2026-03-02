#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>

using namespace dealii;

template <int dim>
struct Geometry
{
  Geometry(const std::string filename);

  Triangulation<dim> triangulation;

  void scale(const double scale_factor);
  void align_line_of_action_with_x_axis();
  void tag_boundaries();
  void write_grid(const std::string filename);
};

template <int dim>
Geometry<dim>::Geometry(const std::string mesh_filename)
{
  GridIn<3> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream infile(mesh_filename);
  if (infile.fail())
    throw std::invalid_argument("Cannot open file: " + mesh_filename +  
          ". Make sure the file exists and it has read permissions.");
  else
    gridin.read_msh(infile);
  infile.close();

  double volume = GridTools::volume(triangulation);
  double density = 1060;
  double mass = density * volume;

  std::cout << "Volume: " << volume << "\n"
            << "Mass:   " << mass << "\n"
            << std::endl;
}

template <int dim>
void Geometry<dim>::scale(const double scale_factor)
{
  GridTools::scale(scale_factor, triangulation);

  double volume = GridTools::volume(triangulation);
  double density = 1060;
  double mass = density * volume;

  std::cout << "Volume (after scaling): " << volume << "\n"
            << "Mass (after scaling):   " << mass << "\n"
            << std::endl;
}

template <int dim>
void Geometry<dim>::align_line_of_action_with_x_axis()
{
  Point<dim> p1(-2.073870e+01, -1.911716e+01, -3.719337e+01);
  Point<dim> p2(-5.061858e+01, 5.845398e+01, 1.950347e+02);
  Tensor<1,dim> line_of_action = p2 - p1;
  
  double theta = std::acos(line_of_action[2]/line_of_action.norm());
  double sign_y = (line_of_action[2] >= 0 ? 1 : -1);
  double phi = sign_y * std::acos(line_of_action[0]/std::sqrt(line_of_action[0] * line_of_action[0] + line_of_action[1] * line_of_action[1]));
  
  GridTools::rotate(Tensor<1,dim>({0,0,1}), -phi, triangulation);
  GridTools::rotate(Tensor<1,dim>({0,1,0}), -theta, triangulation);
  GridTools::rotate(Tensor<1,dim>({0,1,0}), M_PI/2.0, triangulation);
  GridTools::rotate(Tensor<1,dim>({1,0,0}), M_PI, triangulation);
}

template <int dim>
void Geometry<dim>::tag_boundaries()
{
  // Find min_x, max_x coordinates. The difference is the length of the muscle
  typename Triangulation<dim>::vertex_iterator
  vertex = triangulation.begin_vertex(),
  endv = triangulation.end_vertex();

  Point<3> &v0 = vertex->vertex();
  double min_x = v0[0];
  double max_x = v0[0];

  for(; vertex != endv; ++vertex)
  {
    Point<3> &v = vertex->vertex();
    if (v[0] < min_x)
      min_x = v[0];
    if (v[0] > max_x)
      max_x = v[0];
  }

  double muscle_length = max_x - min_x;
  double fraction_left = 0.02;
  double fraction_right = 0.02;
  
  // Loop over all cells and tag boundaries
  for (const auto &face : triangulation.active_face_iterators())
    if (face->at_boundary())
    {
      if (face->center()[0] < min_x + fraction_left * muscle_length)
        face->set_boundary_id(0);
      else if (face->center()[0] > max_x - fraction_right * muscle_length)
        face->set_boundary_id(1);
      else
        face->set_boundary_id(2);
    }
}

template <int dim>
void Geometry<dim>::write_grid(const std::string filename)
{
  GridOut           grid_out;
  GridOutFlags::Msh write_flags;
  write_flags.write_faces = true;
  grid_out.set_flags(write_flags);
  std::ofstream output(filename);
  grid_out.write_msh(triangulation, output);
}

int main()
{
  const unsigned int dim = 3;
  Geometry<dim> mg_mesh("MG_left.msh");
  mg_mesh.align_line_of_action_with_x_axis();
  mg_mesh.scale(0.001);
  mg_mesh.tag_boundaries();
  mg_mesh.write_grid("MG_left_transformed.msh");
  
  return 0;
}