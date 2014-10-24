// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef C_RIDGE_EXTRACTOR_H_
#define C_RIDGE_EXTRACTOR_H_

class vtkStructuredPoints;
class vtkPolyData;

class CRidgeExtractor {
 public:
  static void get_gradient(vtkStructuredPoints *scalar_field,
                           vtkStructuredPoints **gradient);

  static void get_gradient_and_hessian(vtkStructuredPoints *scalar_field,
                                       vtkStructuredPoints **gradient,
                                       vtkStructuredPoints **hessian);

  static void get_cauchy_green_tensor(vtkStructuredPoints *flow_map,
                                      vtkStructuredPoints *cauchy_green);

  // The flow map is stored in a scalar field of flow_map->GetPointData().
  vtkPolyData *extract_ridges(vtkStructuredPoints *flow_map);
};

#endif  // C_RIDGE_EXTRACTOR_H_
