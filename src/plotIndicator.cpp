#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVectorText.h>
#include <vtkVertexGlyphFilter.h>

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

std::vector<std::array<double, 3>> load_points(const std::string &filename) {
  std::vector<std::array<double, 3>> points;
  std::ifstream file(filename);

  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::array<double, 3> point;
      if (iss >> point[0] >> point[1] >> point[2]) {
        points.push_back(point);
      } else {
        std::cerr << "Invalid line in file: " << line << std::endl;
      }
    }
    file.close();
  } else {
    std::cerr << "Error opening file for reading!" << std::endl;
  }

  return points;
}

std::vector<double> loadVectorFromFile(const std::string &filename) {
  std::vector<double> vec;
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return vec;
  }

  std::string line;
  if (std::getline(file, line)) { // Read the single line
    std::stringstream ss(line);
    std::string value;
    while (std::getline(ss, value, ',')) { // Split by comma
      try {
        vec.push_back(std::stod(value)); // Convert to double
      } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid number found: " << value << std::endl;
      }
    }
  }

  file.close();
  return vec;
}

int main(int argc, char **argv) {
  // retrieve values
  std::vector<std::array<double, 3>> vertices = load_points("centers.txt");
  std::vector<std::array<double, 3>> samples = load_points("points.txt");
  std::vector<double> weights = loadVectorFromFile("x_normalized.txt");

  vtkNew<vtkDoubleArray> pointScalars;
  pointScalars->SetName("pointScalars");
  pointScalars->SetNumberOfComponents(1);

  vtkNew<vtkPoints> points;
  for (int i = 0; i < vertices.size(); i++) {
    std::array<double, 3> v = vertices[i];
    points->InsertNextPoint(v[0], v[1], v[2]);
    double scalar = weights[i];
    // double scalar =
    //     sqrt(std::pow(v[0], 2) + std::pow(v[1], 2) + std::pow(v[2], 2)) - 1;
    pointScalars->InsertNextValue(scalar);
    std::cout << scalar << std::endl;
  }

  vtkNew<vtkPoints> points2;
  for (int i = 0; i < samples.size(); i++) {
    std::array<double, 3> v = samples[i];
    points2->InsertNextPoint(v[0], v[1], v[2]);
  }

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
  // vertexActor->GetProperty()->SetSpecularPower(50.0);

  vtkNew<vtkActor> vertexActor2;
  vertexActor2->SetMapper(mapper2);
  vertexActor2->GetProperty()->SetPointSize(3.0);
  vertexActor2->GetProperty()->SetSpecular(1.0);
  // vertexActor2->GetProperty()->SetSpecularPower(50.0);

  // Render
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // int target = points->GetNumberOfPoints() / 10;
  // int i = points->GetNumberOfPoints() - 1;
  // while (target > 0) {
  //   if (weights[i] < 0) {
  //     double p[3];
  //     points->GetPoint(i, p);
  //
  //     vtkNew<vtkVectorText> vectorText;
  //     vectorText->SetText(std::to_string(((int)(weights[i] / 40000)))
  //                             .c_str()); // Use point index as label
  //
  //     vtkNew<vtkTransform> transform;
  //     transform->Translate(p[0], p[1], p[2]);
  //     transform->Scale(0.01, 0.01, 0.01);
  //
  //     vtkNew<vtkTransformPolyDataFilter> transformFilter;
  //     transformFilter->SetTransform(transform);
  //     transformFilter->SetInputConnection(vectorText->GetOutputPort());
  //     transformFilter->Update();
  //
  //     vtkNew<vtkPolyDataMapper> textMapper;
  //     textMapper->SetInputConnection(transformFilter->GetOutputPort());
  //
  //     vtkSmartPointer<vtkActor> textActor = vtkSmartPointer<vtkActor>::New();
  //     textActor->SetMapper(textMapper);
  //
  //     renderer->AddActor(textActor);
  //     target -= 1;
  //   }
  //   i--;
  // }

  // Create a sphere
  vtkNew<vtkSphereSource> sphereSource;
  sphereSource->SetRadius(1.0);         // Unit sphere
  sphereSource->SetThetaResolution(32); // Increase for a smoother sphere
  sphereSource->SetPhiResolution(32);   // Increase for a smoother sphere
  sphereSource->Update();

  // Create a mapper
  vtkNew<vtkPolyDataMapper> sphereMapper;
  sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

  // Create an actor
  vtkNew<vtkActor> sphereActor;
  sphereActor->SetMapper(sphereMapper);
  sphereActor->GetProperty()->SetColor(0.0, 0.0, 0.0);

  renderer->AddActor(vertexActor);
  renderer->AddActor(vertexActor2);
  // renderer->AddActor(sphereActor);

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
