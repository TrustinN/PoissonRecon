#include "utils/plot.hpp"
#include "Normal.hpp"
#include "utils/linalg.hpp"
#include "utils/sampling.hpp"

int main() {
  // std::vector<std::array<double, 3>> vertices = sample_sphere(10000, 3);
  std::vector<std::array<double, 3>> vertices = sample_box(10000, 1.5, 2, 1);
  NormalApproximations na(vertices);
  std::vector<std::array<double, 3>> normals = na.normals();

  vtkNew<vtkPoints> points = load_points(vertices);

  std::vector<std::array<double, 3>> b(vertices.size());
  std::vector<double> weights(vertices.size());
  for (int i = 0; i < vertices.size(); i++) {
    b[i] = vertices[i] + normals[i];
    weights[i] = (dot(normals[i], {0.3, 0.3, 0.3}) + 1) / 2;
  }
  vtkNew<vtkCellArray> linesArray = load_lines(vertices, b, points);
  vtkNew<vtkDoubleArray> lineScalars = load_scalars(weights);

  // data here
  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(points);
  polyData->SetLines(linesArray);
  polyData->GetCellData()->SetScalars(lineScalars);

  vtkNew<vtkLookupTable> rainbowBlueRedLut;
  rainbowBlueRedLut->SetNumberOfColors(256);
  rainbowBlueRedLut->SetHueRange(0.667, 0.0);
  rainbowBlueRedLut->Build();

  // Mapper and Actor for Lines
  vtkNew<vtkPolyDataMapper> lineMapper;
  lineMapper->SetInputData(polyData);
  lineMapper->SetLookupTable(rainbowBlueRedLut);
  lineMapper->SetScalarRange(lineScalars->GetRange());

  vtkNew<vtkActor> lineActor;
  lineActor->SetMapper(lineMapper);
  lineActor->GetProperty()->SetLineWidth(0.01);
  lineActor->GetProperty()->SetInterpolationToPhong();

  // vertex Actor & Mapper
  vtkNew<vtkVertexGlyphFilter> vertexFilter;
  vertexFilter->SetInputData(polyData);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(vertexFilter->GetOutputPort());

  vtkNew<vtkActor> vertexActor;
  vertexActor->SetMapper(mapper);
  vertexActor->GetProperty()->SetPointSize(3.3);
  vertexActor->GetProperty()->SetSpecular(1.0);
  vertexActor->GetProperty()->SetSpecularPower(50.0);

  // Render
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(vertexActor);
  renderer->AddActor(lineActor);

  // Window Render
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}
