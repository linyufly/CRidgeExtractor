// Author: Mingcheng Chen (linyufly@gmail.com)

#include "c_ridge_extractor.h"
#include "util.h"
#include "math.h"
#include "marchingCubesTable.h"

#include <vtkStructuredPoints.h>
#include <vtkPolyData.h>
#include <vtkGradientFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkMath.h>

#include <cstring>

#include <algorithm>

namespace {

const int kVertexList[8][3] = {
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
};

const int kEdgeList[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
    {4, 5}, {5, 6}, {6, 7}, {7, 4},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

double dot_product_3d(double *a, double *b) {
  double result = 0.0;
  for (int i = 0; i < 3; i++) {
    result += a[i] * b[i];
  }
  return result;
}

// e3 is the eigen vector of the smallest eigen value.
// l3 is the smallest eigen value.
void get_e3(double **hessian, double *e3, double *l3) {
  double *eigen_values = new double[3];
  double **eigen_vectors = create_matrix<double>(3, 3);
  vtkMath::Jacobi(hessian, eigen_values, eigen_vectors);

  *l3 = eigen_values[2];

  for (int i = 0; i < 3; i++) {
    e3[i] = eigen_vectors[i][2];
  }

  delete [] eigen_values;
  delete_matrix(eigen_vectors);
  delete_matrix(hessian);
}

}

void CRidgeExtractor::get_gradient(
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints **gradient) {
  vtkSmartPointer<vtkGradientFilter> grad_filter =
      vtkSmartPointer<vtkGradientFilter>::New();
  grad_filter->SetInputData(scalar_field);
  grad_filter->Update();

  *gradient = vtkStructuredPoints::New();
  (*gradient)->ShallowCopy(grad_filter->GetOutput());

  (*gradient)->GetPointData()->SetActiveScalars("Gradients");
  (*gradient)->GetPointData()->GetScalars()->SetName("gradient");
}

void CRidgeExtractor::get_gradient_and_hessian(
    vtkStructuredPoints *scalar_field,
    vtkStructuredPoints **gradient,
    vtkStructuredPoints **hessian) {
  // Calculate the gradient
  vtkSmartPointer<vtkGradientFilter> grad_filter =
      vtkSmartPointer<vtkGradientFilter>::New();
  grad_filter->SetInputData(scalar_field);
  grad_filter->Update();

  *gradient = vtkStructuredPoints::New();
  (*gradient)->ShallowCopy(grad_filter->GetOutput());

  (*gradient)->GetPointData()->SetActiveScalars("Gradients");
  (*gradient)->GetPointData()->GetScalars()->SetName("gradient");

  // Calculate the hessian
  grad_filter = vtkSmartPointer<vtkGradientFilter>::New();
  grad_filter->SetInputData(*gradient);
  grad_filter->Update();

  *hessian = vtkStructuredPoints::New();
  (*hessian)->ShallowCopy(grad_filter->GetOutput());

  (*hessian)->GetPointData()->SetActiveScalars("Gradients");
  (*hessian)->GetPointData()->GetScalars()->SetName("hessian");
}

void CRidgeExtractor::get_cauchy_green_tensor(
    vtkStructuredPoints *flow_map, vtkStructuredPoints **cauchy_green) {
  get_gradient(flow_map, cauchy_green);
  (*cauchy_green)->GetPointData()->GetScalars()->SetName("cauchy_green");
}

vtkPolyData *CRidgeExtractor::extract_ridges(
    vtkStructuredPoints *flow_map) {
  int dimensions[3];
  double spacing[3], origin[3];
  scalar_field->GetDimensions(dimensions);
  scalar_field->GetSpacing(spacing);
  scalar_field->GetOrigin(origin);

  // Calculate the gradient and hessian
  vtkStructuredPoints *gradient_field = NULL, *hessian_field = NULL;
  get_gradient_and_hessian(scalar_field, &gradient_field, &hessian_field);

  int nx = dimensions[0];
  int ny = dimensions[1];
  int nz = dimensions[2];

  int ****edge_mark = create_4d_array<int>(nx, ny, nz, 3);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
        for (int d = 0; d < 3; d++) {
          edge_mark[x][y][z][d] = -1;
        }
      }
    }
  }

  vtkSmartPointer<vtkPoints> mesh_points =
      vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> mesh_cells =
      vtkSmartPointer<vtkCellArray>::New();

  for (int x = 0; x + 1 < nx; x++) {
    for (int y = 0; y + 1 < ny; y++) {
      for (int z = 0; z + 1 < nz; z++) {
        double dot_prod[3][3][3], e3[3][3][3][3], grad[3][3][3][3];

        // Collect e3 and grad
        int num_pos = 0;

        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              int curr_x = x + dx;
              int curr_y = y + dy;
              int curr_z = z + dz;

              double **hessian = create_matrix<double>(3, 3);
              int point_id = (curr_z * ny + curr_y) * nx + curr_x;

              double tensor[9];
              hessian_field->GetPointData()->GetScalars()
                                           ->GetTuple(point_id, tensor);
              for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                  hessian[i][j] = tensor[i * 3 + j];
                }
              }

              double l3;
              get_e3(hessian, e3[dx][dy][dz], &l3);
              if (l3 > 0.0) {
                num_pos++;
              }
              
              gradient_field->GetPointData()->GetScalars()
                                            ->GetTuple(point_id, tensor);
              for (int i = 0; i < 3; i++) {
                grad[dx][dy][dz][i] = tensor[i];
              }
            }
          }
        }

        if (num_pos > 0) {
          continue;
        }

        // Re-orientate e3
        double **vectors = create_matrix<double>(8, 3);
        int num_vectors = 0;
        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              for (int i = 0; i < 3; i++) {
                vectors[num_vectors][i] = e3[dx][dy][dz][i];
              }
              num_vectors++;
            }
          }
        }

        double *pivot = principal_component(vectors, 8, 3);

        delete_matrix(vectors);

        for (int dx = 0; dx < 2; dx++) {
          for (int dy = 0; dy < 2; dy++) {
            for (int dz = 0; dz < 2; dz++) {
              if (dot_product_3d(pivot, e3[dx][dy][dz]) < 0.0) {
                for (int i = 0; i < 3; i++) {
                  e3[dx][dy][dz][i] *= -1.0;
                }
              }

              dot_prod[dx][dy][dz] = dot_product_3d(grad[dx][dy][dz],
                                                    e3[dx][dy][dz]);
            }
          }
        }

        delete [] pivot;

        // Identify iso-surfaces
        int cube_code = 0;
        for (int i = 0; i < 8; i++) {
          int dx = kVertexList[i][0];
          int dy = kVertexList[i][1];
          int dz = kVertexList[i][2];

          if (dot_prod[dx][dy][dz] <= 0.0) {
            cube_code |= (1 << i);
          }
        }

        for (int i = 0; i < numVertsTable[cube_code]; i += 3) {
          mesh_cells->InsertNextCell(3);

          for (int j = 0; j < 3; j++) {
            int edge_idx = triTable[cube_code][i + j];
            int vtx_1 = kEdgeList[edge_idx][0];
            int vtx_2 = kEdgeList[edge_idx][1];

            int dim;
            for (dim = 0; dim < 3; dim++) {
              if (kVertexList[vtx_1][dim] != kVertexList[vtx_2][dim]) {
                break;
              }
            }

            if (kVertexList[vtx_1][dim] > kVertexList[vtx_2][dim]) {
              std::swap(vtx_1, vtx_2);
            }

            int start_x = kVertexList[vtx_1][0];
            int start_y = kVertexList[vtx_1][1];
            int start_z = kVertexList[vtx_1][2];

            int finish_x = kVertexList[vtx_2][0];
            int finish_y = kVertexList[vtx_2][1];
            int finish_z = kVertexList[vtx_2][2];

            // Insert a new point to the mesh if necessary
            if (edge_mark[x + start_x][y + start_y][z + start_z][dim] == -1) {
              edge_mark[x + start_x][y + start_y][z + start_z][dim] =
                  mesh_points->GetNumberOfPoints();

              double dot_prod_1 = dot_prod[start_x][start_y][start_z];
              double dot_prod_2 = dot_prod[finish_x][finish_y][finish_z];

              if (dot_prod_1 * dot_prod_2 > 0.0) {
                report_error("Same sign in marching cubes");
              }

              double lambda = -dot_prod_1 / (dot_prod_2 - dot_prod_1);

              double aug_x = (finish_x - start_x) * spacing[0] * lambda;
              double aug_y = (finish_y - start_y) * spacing[1] * lambda;
              double aug_z = (finish_z - start_z) * spacing[2] * lambda;

              double point_x = origin[0] + spacing[0] * (x + start_x) + aug_x;
              double point_y = origin[1] + spacing[1] * (y + start_y) + aug_y;
              double point_z = origin[2] + spacing[2] * (z + start_z) + aug_z;

              mesh_points->InsertNextPoint(point_x, point_y, point_z);
            }

            mesh_cells->InsertCellPoint(
                edge_mark[x + start_x][y + start_y][z + start_z][dim]);
          }
        }
      }
    }
  }

  vtkPolyData *mesh = vtkPolyData::New();
  mesh->SetPoints(mesh_points);
  mesh->SetPolys(mesh_cells);

  delete_4d_array(edge_mark);

  gradient_field->Delete();
  hessian_field->Delete();

  return mesh;
}
