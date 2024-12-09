#include "BSpline.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>

std::array<double, 6> generateBounds() {
  std::array<double, 6> bounds = {-1.5, 1.5, -1.5, 1.5, -1.5, 1.5};
  return bounds;
}

std::array<double, 6> generateBounds(const std::array<double, 3> &center,
                                     double width) {
  std::array<double, 6> bounds = {-1.5, 1.5, -1.5, 1.5, -1.5, 1.5};
  std::array<double, 6> shift = {center[0], center[0], center[1],
                                 center[1], center[2], center[2]};
  bounds = width * bounds + shift;
  return bounds;
}

struct sfData {
  std::array<double, 6> bounds;
  ScalarField<2> sf;

  sfData() {
    bounds = generateBounds();
    sf = ScalarField<2>(BSpline);
  }

  sfData(double width) {
    std::array<double, 3> c = {0, 0, 0};
    bounds = generateBounds(c, width);
    sf = ScalarField<2>(BSpline, c, 1.0 / width);
  }

  sfData(std::array<double, 3> center) {
    bounds = generateBounds(center, 1.0);
    sf = ScalarField<2>(BSpline, center);
  }

  sfData(std::array<double, 3> center, double w) {
    bounds = generateBounds(center, w);
    sf = ScalarField<2>(BSpline, center, 1.0 / w);
  }
};

void renderScalarField(sfData sf, double weight, double factor,
                       int resolution[3],
                       vtkSmartPointer<vtkRenderer> renderer) {
  auto bounds = sf.bounds;
  auto scalarf = sf.sf;

  // Create the data structure
  vtkSmartPointer<vtkStructuredPoints> structuredPoints =
      vtkSmartPointer<vtkStructuredPoints>::New();
  structuredPoints->SetDimensions(resolution[0] - 1, resolution[1] - 1,
                                  resolution[2] - 1);
  structuredPoints->SetOrigin(bounds[0], bounds[2], bounds[4]);
  structuredPoints->SetSpacing((bounds[1] - bounds[0]) / (resolution[0] - 1),
                               (bounds[3] - bounds[2]) / (resolution[1] - 1),
                               (bounds[5] - bounds[4]) / (resolution[2] - 1));

  // Populate the scalar field
  vtkSmartPointer<vtkDoubleArray> scalars =
      vtkSmartPointer<vtkDoubleArray>::New();
  int dims[3];
  structuredPoints->GetDimensions(dims);

  for (int z = 0; z < dims[2]; ++z) {
    for (int y = 0; y < dims[1]; ++y) {
      for (int x = 0; x < dims[0]; ++x) {
        // Calculate the physical coordinates based on the bounds and resolution
        double coordX = bounds[0] + x * (bounds[1] - bounds[0]) / (dims[0]);
        double coordY = bounds[2] + y * (bounds[3] - bounds[2]) / (dims[1]);
        double coordZ = bounds[4] + z * (bounds[5] - bounds[4]) / (dims[2]);

        // Evaluate the scalar field at the calculated coordinates
        double scalar = factor * weight * scalarf({coordX, coordY, coordZ});

        // Insert the scalar value into the array
        scalars->InsertNextValue(scalar);
      }
    }
  }
  structuredPoints->GetPointData()->SetScalars(scalars);

  // Create a color transfer function to map scalar values to colors
  vtkSmartPointer<vtkColorTransferFunction> colorFunc =
      vtkSmartPointer<vtkColorTransferFunction>::New();

  // Add multiple control points for smoother transitions
  colorFunc->AddRGBPoint(-1.0, 0.0, 0.0, 1.0); // Blue for minimum value
  colorFunc->AddRGBPoint(-0.5, 0.0, 0.5, 1.0); // Light blue
  colorFunc->AddRGBPoint(0.0, 0.0, 1.0, 0.0);  // Green for zero
  colorFunc->AddRGBPoint(0.5, 1.0, 0.5, 0.0);  // Orange
  colorFunc->AddRGBPoint(1.0, 1.0, 0.0, 0.0);  // Red for maximum value

  // Generate contours
  vtkSmartPointer<vtkContourFilter> contourFilter =
      vtkSmartPointer<vtkContourFilter>::New();
  contourFilter->SetInputData(structuredPoints);
  contourFilter->SetNumberOfContours(100);
  vtkSmartPointer<vtkDataArray> dataArray = scalars;
  double range[2];
  dataArray->GetRange(range); // Get the scalar range
  double step =
      (range[1] - range[0]) / 100; // Divide the range into 100 contours

  for (int i = 0; i < 100; ++i) {
    contourFilter->SetValue(
        i, range[0] + i * step); // Set contour levels based on the scalar range
  }

  // Set up the mapper
  vtkSmartPointer<vtkDataSetMapper> mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputConnection(contourFilter->GetOutputPort());
  mapper->SetLookupTable(colorFunc);

  // Set up the actor
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetOpacity(0.1);

  renderer->AddActor(actor);
}

int main(int argc, char **argv) {
  int depth = 4;
  int end = -1;
  bool iso = false;
  double factor = 1.0;
  if (argc > 1) {
    for (int i = 0; i < argc; i++) {
      std::string option = std::string(argv[i]);
      if (option == "--single") {
        i += 1;
        try {
          depth = std::stoi(argv[i]);
          end = depth - 1;
        } catch (std::exception) {
          i -= 1;
        }
      } else if (option == "--cascade") {
        try {
          i += 1;
          depth = std::stoi(argv[i]);
          i += 1;
          end = std::stoi(argv[i]);
        } catch (std::exception) {
          std::cout << "error" << std::endl;
        }
      } else if (option == "--iso") {
        iso = true;
      } else if (option == "--factor") {
        i += 1;
        factor = std::stof(argv[i]);
        std::cout << factor << std::endl;
      }
    }
  }
  int resolution1[3] = {10, 10, 10};

  // Set up the renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

  for (int i = depth; i > end; i--) {
    // retrieve values
    std::cout << "Starting level: " << i << std::endl;

    std::vector<std::array<double, 3>> vertices =
        load_points("data/centers_depth_" + std::to_string(i) + ".txt");
    std::vector<double> weights;
    if (!iso) {
      weights =
          loadVectorFromFile("data/x_depth_" + std::to_string(i) + ".txt");
    } else {
      weights = loadVectorFromFile("data/iso_vals_depth_" + std::to_string(i) +
                                   ".txt");
    }
    std::vector<double> widths =
        loadVectorFromFile("data/widths_depth_" + std::to_string(i) + ".txt");

    for (int j = 0; j < vertices.size(); j++) {
      sfData sf(vertices[j], widths[j]);
      renderScalarField(sf, weights[j], factor, resolution1, renderer);
    }
  }

  // Set up the render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // Set up the interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Start the rendering loop
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}
