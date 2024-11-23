#include "Normal.hpp"
#include "utils/linalg.hpp"
#include "utils/sampling.hpp"
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkVertexGlyphFilter.h>

int main() {
  // std::vector<std::array<double, 3>> vertices = sample_sphere(10000, 3);
  std::vector<std::array<double, 3>> vertices = sample_box(10000, 1.5, 2, 1);
  NormalApproximations na(vertices);
  std::vector<std::array<double, 3>> normals = na.normals();

  vtkNew<vtkPoints> points;
  for (const auto &v : vertices) {
    points->InsertNextPoint(v[0], v[1], v[2]);
  }

  vtkNew<vtkCellArray> linesArray;
  vtkNew<vtkDoubleArray> lineScalars;
  lineScalars->SetName("LineScalars");
  lineScalars->SetNumberOfComponents(1);

  // Add lines to points and cell array
  for (int i = 0; i < normals.size(); i++) {
    std::array<double, 3> n = normals[i];
    std::array<double, 3> extend = {
        vertices[i][0] + n[0], vertices[i][1] + n[1], vertices[i][2] + n[2]};
    points->InsertNextPoint(extend.data());
    points->InsertNextPoint(vertices[i].data());

    vtkNew<vtkLine> vtkLine;
    vtkLine->GetPointIds()->SetId(0, points->GetNumberOfPoints() - 2);
    vtkLine->GetPointIds()->SetId(1, points->GetNumberOfPoints() - 1);
    linesArray->InsertNextCell(vtkLine);

    // Assign a scalar value to the line
    double scalar = (dot(n, {0.3, 0.3, 0.3}) + 1) / 2;
    lineScalars->InsertNextValue(scalar);
  }

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

  // Get the bounds of the points (minX, maxX, minY, maxY, minZ, maxZ)
  double bounds[6];
  polyData->GetBounds(bounds);

  // Calculate the center of the points
  double centerX = (bounds[0] + bounds[1]) / 2.0;
  double centerY = (bounds[2] + bounds[3]) / 2.0;
  double centerZ = (bounds[4] + bounds[5]) / 2.0;

  vtkSmartPointer<vtkCamera> camera = renderer->GetActiveCamera();
  camera->SetPosition(centerX, centerY, bounds[5] + 5);
  camera->SetFocalPoint(centerX, centerY, centerZ);
  camera->SetViewUp(0, 1, 0);

  // Window Render
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}
