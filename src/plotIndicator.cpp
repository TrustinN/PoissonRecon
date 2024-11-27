#include "utils/io.hpp"
#include "utils/plot.hpp"
#include <array>
#include <vector>

int main(int argc, char **argv) {
  // retrieve values
  std::vector<std::array<double, 3>> vertices = load_points("centers.txt");
  std::vector<std::array<double, 3>> samples = load_points("points.txt");
  std::vector<double> weights = loadVectorFromFile("x_normalized.txt");

  vtkNew<vtkPoints> points = load_points(vertices);
  vtkNew<vtkPoints> points2 = load_points(samples);
  vtkNew<vtkDoubleArray> pointScalars = load_scalars(weights);

  double range[2];
  pointScalars->GetRange(range);
  std::cout << "Scalar Range: [" << range[0] << ", " << range[1] << "]"
            << std::endl;

  // data here
  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(points);
  polyData->GetPointData()->SetScalars(pointScalars);

  vtkNew<vtkPolyData> polyData2;
  polyData2->SetPoints(points2);

  vtkNew<vtkLookupTable> rainbowBlueRedLut;
  rainbowBlueRedLut->SetNumberOfColors(256);
  rainbowBlueRedLut->SetHueRange(0.667, 0);
  rainbowBlueRedLut->Build();

  // vertex Actor & Mapper
  vtkNew<vtkVertexGlyphFilter> vertexFilter;
  vertexFilter->SetInputData(polyData);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(vertexFilter->GetOutputPort());
  mapper->SetLookupTable(rainbowBlueRedLut);
  mapper->SetScalarRange(-.1, .1);
  mapper->SetColorModeToMapScalars();
  mapper->SetScalarModeToUsePointData();

  vtkNew<vtkVertexGlyphFilter> vertexFilter2;
  vertexFilter2->SetInputData(polyData2);

  vtkNew<vtkPolyDataMapper> mapper2;
  mapper2->SetInputConnection(vertexFilter2->GetOutputPort());

  vtkNew<vtkActor> vertexActor;
  vertexActor->SetMapper(mapper);
  vertexActor->GetProperty()->SetPointSize(3.0);
  vertexActor->GetProperty()->SetSpecular(1.0);

  vtkNew<vtkActor> vertexActor2;
  vertexActor2->SetMapper(mapper2);
  vertexActor2->GetProperty()->SetPointSize(3.0);
  vertexActor2->GetProperty()->SetSpecular(1.0);

  // Render
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(vertexActor);
  renderer->AddActor(vertexActor2);

  // Window Render
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}
