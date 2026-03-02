/*
 * ============================================================================
 *  Subject-specific Mesh Preprocessor (deal.II-based)
 * ============================================================================
 *
 *  IMPORTANT
 *  ---------
 *  This preprocessor is subject-specific and is designed to work ONLY
 *  with the mesh file:
 *
 *      MG_left.msh
 *
 *  The geometric alignment parameters and anatomical reference points are
 *  hard-coded for this particular mesh. Using a different mesh without
 *  modifying the code will produce incorrect alignment and boundary tags.
 *
 *
 *  Description
 *  -----------
 *  This program reads the MG_left.msh 3D Gmsh mesh, applies a sequence of
 *  geometric transformations, assigns boundary IDs based on anatomical
 *  position, and writes a transformed mesh:
 *
 *      MG_left_transformed.msh
 *
 *  The tool is intended for preprocessing this specific skeletal muscle
 *  geometry prior to Flexodeal simulations.
 *
 *
 *  Workflow
 *  --------
 *  The program performs the following steps:
 *
 *    1. Load mesh from MG_left.msh.
 *    2. Compute and print:
 *         - Volume
 *         - Estimated mass (density = 1060 kg/m^3).
 *    3. Align a predefined line of action with the global x-axis.
 *    4. Apply uniform scaling (currently hard-coded to 0.001).
 *    5. Automatically tag boundary faces.
 *    6. Write the transformed mesh to MG_left_transformed.msh.
 *
 *
 *  Boundary Tagging Convention
 *  ----------------------------
 *  After alignment, the muscle is assumed to extend primarily along the
 *  x-direction. Boundary IDs are assigned as follows:
 *
 *      boundary_id = 0  →  Left end (first 2% of muscle length)
 *      boundary_id = 1  →  Right end (last 2% of muscle length)
 *      boundary_id = 2  →  Lateral surface (remaining boundary faces)
 *
 *  These IDs are intended for later application of boundary conditions in
 *  finite element simulations.
 *
 *
 *  Alignment Procedure
 *  --------------------
 *  The function align_line_of_action_with_x_axis() rotates the mesh so that
 *  a predefined line (between two hard-coded anatomical points) becomes
 *  aligned with the global x-axis.
 *
 *  Because these points are specific to MG_left.msh, this transformation
 *  is NOT general-purpose.
 *
 *
 *  Units
 *  -----
 *  - The original mesh is assumed to be in millimetres.
 *  - A uniform scale factor of 0.001 is applied (mm → m).
 *  - Density is assumed to be 1060 kg/m^3 (typical skeletal muscle value).
 *
 *
 *  Limitations
 *  -----------
 *  - Works ONLY for MG_left.msh.
 *  - Not intended as a general mesh preprocessor.
 *  - Alignment parameters are hard-coded.
 *  - Boundary tagging assumes the principal axis is x after alignment.
 *  - Designed for 3D meshes (dim = 3).
 *
 *
 *  Author
 *  ------------
 *  Javier Almonacid
 *  Research Assistant
 *  Neuromuscular Mechanics Laboratory
 *
 * ============================================================================
 */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>

using namespace dealii;

// We first create a simple class to read, store, transform, and output the
// triangulation (mesh).
template <int dim>
struct Geometry
{
  Geometry(const std::string filename);

  Triangulation<dim> triangulation;

  // The chosen mesh requires these 3 operations
  void scale(const double scale_factor);
  void align_line_of_action_with_x_axis();
  void tag_boundaries();
  // A simple function to output the transformed triangulation.
  void write_grid(const std::string filename);
};

// The constructor is in charge of reading the triangulation
// and computing the volume and mass of the muscle.
// The result is given in the original mesh's length units.
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

  std::cout << "\nMuscle density:          " << density << " kg/m^3\n";

  std::cout << "\nVolume (before scaling): " << volume << "\n"
            << "Mass (before scaling):   " << mass << "\n"
            << std::endl;
}

// This function scales the triangulation. For example, if the mesh
// is given in milimetres, then we apply a scale factor of 0.001 to
// obtain a mesh in metres and a mass in kg.
template <int dim>
void Geometry<dim>::scale(const double scale_factor)
{
  GridTools::scale(scale_factor, triangulation);

  double volume = GridTools::volume(triangulation);
  double density = 1060;
  double mass = density * volume;

  std::cout << "Volume (after scaling):  " << volume << "\n"
            << "Mass (after scaling):    " << mass << "\n"
            << std::endl;
}

// To construct a line of action, we choose two mesh vertices, one from
// each end of the muscle. The difference of these 2 points creates the 
// vector line_of_action. Note that this choice of two points is not
// unique.
template <int dim>
void Geometry<dim>::align_line_of_action_with_x_axis()
{
  Point<dim> p1(-2.073870e+01, -1.911716e+01, -3.719337e+01);
  Point<dim> p2(-5.061858e+01, 5.845398e+01, 1.950347e+02);
  Tensor<1,dim> line_of_action = p2 - p1;
  
  // -------------------------------------------------------------------------
  // STEP 1:
  // Compute polar angle (theta) between the line of action and the z-axis.
  //
  //   cos(theta) = v_z / ||v||
  //
  // After rotating by -theta about the y-axis, the vector will lie in the
  // x–z plane (i.e., remove its inclination relative to z).
  // -------------------------------------------------------------------------
  double theta = std::acos(line_of_action[2]/line_of_action.norm());
  // Determine sign for azimuthal angle to reduce quadrant ambiguity
  double sign_y = (line_of_action[2] >= 0 ? 1 : -1);
  // -------------------------------------------------------------------------
  // STEP 2:
  // Compute azimuthal angle (phi) of the projection of the vector onto
  // the x–y plane.
  //
  //   cos(phi) = v_x / sqrt(v_x^2 + v_y^2)
  //
  // After rotating by -phi about the z-axis, the projection of the vector
  // onto the x–y plane aligns with the x-axis.
  // -------------------------------------------------------------------------
  double phi = sign_y * 
               std::acos(line_of_action[0] / 
               std::sqrt(line_of_action[0] * line_of_action[0] + line_of_action[1] * line_of_action[1]));

  // -------------------------------------------------------------------------
  // Apply rotations to the entire triangulation
  // -------------------------------------------------------------------------
  
  // 1) Remove azimuthal angle (rotate about global z-axis)
  GridTools::rotate(Tensor<1,dim>({0,0,1}), -phi, triangulation);
  // 2) Remove polar inclination (rotate about global y-axis)
  GridTools::rotate(Tensor<1,dim>({0,1,0}), -theta, triangulation);
  // 3) Rotate 90 degrees about y-axis to map aligned z-direction to x-axis
  GridTools::rotate(Tensor<1,dim>({0,1,0}), M_PI/2.0, triangulation);
  // 4) Rotate 180 degrees about x-axis to enforce consistent orientation
  GridTools::rotate(Tensor<1,dim>({1,0,0}), M_PI, triangulation);
}

// This function tags the boundaries *assuming the mesh has been aligned with
// the x-axis*.
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
  // Decide the fraction of the domain (in the longitudinal direction)
  // that will be considered a Dirichlet boundary.
  double fraction_left = 0.02;
  double fraction_right = 0.02;
  
  // Loop over all cells and tag boundaries. We want:
  // - Left end with boundary ID = 0
  // - Right end with boundary ID = 1
  // - All other faces with boundary ID = 2
  // This way, we can match the default setting in Flexodeal:
  // - Boundary ID 0 will be fixed (no displacement),
  // - Boundary ID 1 will have a prescribed strain,
  // - Boundary ID 2 will be interpreted as a traction-free surface.
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

// This function simply writes the transformed mesh with the 
// provided filename.
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

// Main runner. Note that for simplicity, all arguments are hard-coded.
// This means that every time you modify these arguments, you must "make"
// this file again. Ideally, the arguments below should be parsed as
// command-line arguments, but for this tutorial-style file, hard-coding
// them should suffice.
int main()
{
  const unsigned int dim = 3;
  Geometry<dim> mg_mesh("../MG_left.msh");
  mg_mesh.align_line_of_action_with_x_axis();
  mg_mesh.scale(0.001);
  mg_mesh.tag_boundaries();
  mg_mesh.write_grid("../MG_left_transformed.msh");
  
  return 0;
}