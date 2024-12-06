#include "utils/io.hpp"
#include "utils/plot.hpp"
#include <array>
#include <vector>

int main(int argc, char **argv) {

  int depth = 6;
  int end = -1;
  if (argc > 1) {
    depth = std::stoi(argv[1]);
    end = depth - 1;
  }

  std::vector<std::array<double, 3>> samples = load_points("points.txt");
  vtkNew<vtkPoints> points2 = load_points(samples);

  vtkNew<vtkPoints> points;
  vtkNew<vtkDoubleArray> pointScalars;
  pointScalars->SetNumberOfComponents(1);

  for (int i = depth; i > end; i--) {
    // retrieve values

    std::vector<std::array<double, 3>> vertices =
        load_points("data/centers_depth_" + std::to_string(i) + ".txt");
    std::vector<double> weights =
        loadVectorFromFile("data/x_depth_" + std::to_string(i) + ".txt");

    for (int j = 0; j < vertices.size(); j++) {
      std::array<double, 3> p = vertices[j];
      points->InsertNextPoint(p[0], p[1], p[2]);
    }

    for (int j = 0; j < weights.size(); j++) {
      pointScalars->InsertNextValue(weights[j]);
    }

    double range[2];
    pointScalars->GetRange(range);
    std::cout << "Scalar Range: [" << range[0] << ", " << range[1] << "]"
              << std::endl;
  }

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
  mapper->SetScalarRange(-.3, .3);
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
