#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include "utils/plot.hpp"
#include "utils/sampling.hpp"

int main() {
  std::vector<std::array<double, 3>> vertices = load_points("debug/points.txt");
  std::vector<std::array<double, 3>> normals = load_points("debug/normals.txt");
  vtkNew<vtkPoints> points = load_points(vertices);

  std::vector<std::array<double, 3>> b(vertices.size());
  std::vector<double> weights = loadVectorFromFile("debug/weights.txt");
  for (int i = 0; i < vertices.size(); i++) {
    weights[i] *= 10000;
    b[i] = vertices[i] + normals[i] * 4;
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

  std::vector<std::array<double, 3>> v2 = load_points("debug/centers.txt");
  std::vector<double> coeff = loadVectorFromFile("debug/coeff.txt");
  vtkNew<vtkDoubleArray> coeffScalars = load_scalars(coeff);
  vtkNew<vtkPoints> p2 = load_points(v2);
  vtkNew<vtkPolyData> pd2;
  pd2->SetPoints(p2);
  pd2->GetPointData()->SetScalars(coeffScalars);
  vtkNew<vtkVertexGlyphFilter> vf;
  vf->SetInputData(pd2);
  vtkNew<vtkPolyDataMapper> mp2;
  mp2->SetInputConnection(vf->GetOutputPort());
  mp2->SetLookupTable(rainbowBlueRedLut);
  vtkNew<vtkActor> va;
  va->SetMapper(mp2);
  va->GetProperty()->SetPointSize(10);
  va->GetProperty()->SetSpecular(1.0);
  va->GetProperty()->SetSpecularPower(50.0);

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
  renderer->AddActor(va);
  renderer->AddActor(lineActor);

  // Window Render
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}
