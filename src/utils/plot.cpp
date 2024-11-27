#include "plot.hpp"

vtkNew<vtkPoints>
load_points(const std::vector<std::array<double, 3>> &points) {
  vtkNew<vtkPoints> v_points;
  for (const auto &p : points) {
    v_points->InsertNextPoint(p[0], p[1], p[2]);
  }
  return v_points;
};

vtkNew<vtkDoubleArray> load_scalars(const std::vector<double> &weights) {
  vtkNew<vtkDoubleArray> scalars;
  scalars->SetNumberOfComponents(1);
  for (const auto &w : weights) {
    scalars->InsertNextValue(w);
  }
  return scalars;
};

vtkNew<vtkCellArray> load_lines(const std::vector<std::array<double, 3>> &a,
                                const std::vector<std::array<double, 3>> &b,
                                vtkNew<vtkPoints> &points) {
  vtkNew<vtkCellArray> linesArray;

  for (int i = 0; i < a.size(); i++) {
    points->InsertNextPoint(a[i][0], a[i][1], a[i][2]);
    points->InsertNextPoint(b[i][0], b[i][1], b[i][2]);

    vtkNew<vtkLine> line;
    line->GetPointIds()->SetId(0, points->GetNumberOfPoints() - 2);
    line->GetPointIds()->SetId(1, points->GetNumberOfPoints() - 1);
    linesArray->InsertNextCell(line);
  }

  return linesArray;
};
