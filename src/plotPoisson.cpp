#include "Normal.hpp"
#include "PoissonRecon.hpp"
#include "utils/linalg.hpp"
#include "utils/plot.hpp"
#include "utils/sampling.hpp"

void render_points(vtkNew<vtkPoints> &points, vtkNew<vtkRenderer> &renderer) {
  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(points);

  // vertex Actor & Mapper
  vtkNew<vtkVertexGlyphFilter> vertexFilter;
  vertexFilter->SetInputData(polyData);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(vertexFilter->GetOutputPort());

  vtkNew<vtkActor> vertexActor;
  vertexActor->SetMapper(mapper);
  vertexActor->GetProperty()->SetPointSize(3.5);
  vertexActor->GetProperty()->SetSpecular(1.0);
  vertexActor->GetProperty()->SetSpecularPower(50.0);

  renderer->AddActor(vertexActor);
}

void render_lines(vtkNew<vtkCellArray> &lines, vtkNew<vtkDoubleArray> &scalars,
                  vtkNew<vtkPoints> &points, vtkNew<vtkRenderer> &renderer) {
  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(points);
  polyData->SetLines(lines);
  polyData->GetCellData()->SetScalars(scalars);

  vtkNew<vtkLookupTable> rainbowBlueRedLut;
  rainbowBlueRedLut->SetNumberOfColors(256);
  rainbowBlueRedLut->SetHueRange(0.667, 0.0);
  rainbowBlueRedLut->Build();

  // Mapper and Actor for Lines
  vtkNew<vtkPolyDataMapper> lineMapper;
  lineMapper->SetInputData(polyData);
  lineMapper->SetLookupTable(rainbowBlueRedLut);
  lineMapper->SetScalarRange(scalars->GetRange());

  vtkNew<vtkActor> lineActor;
  lineActor->SetMapper(lineMapper);
  lineActor->GetProperty()->SetLineWidth(0.01);
  lineActor->GetProperty()->SetInterpolationToPhong();

  renderer->AddActor(lineActor);
}

int main(int argc, char **argv) {
  std::string obj = argv[1];
  std::vector<std::array<double, 3>> vertices;
  if (obj == "ball") {
    vertices = sample_sphere(1000000, 3);
  } else if (obj == "box") {
    vertices = sample_box(1000000, 1.5, 2, 1);
  }

  NormalApproximations na(vertices);
  PoissonRecon poisson(vertices, na.normals(), na.inward_normals(), 6);
  std::cout << "OK" << std::endl;

  // auto field_nodes = poisson.octree().getNodesAtDepth(6);
  // std::cout << "Getting Nodes" << std::endl;
  // std::vector<std::array<double, 3>> field_points;
  // for (Node *node : field_nodes) {
  //   assert(node != nullptr);
  //   field_points.push_back(node->center);
  // }

  auto field_vertices = poisson.field_centers();
  auto field_normals = poisson.field_normals();

  // // plotting here
  vtkNew<vtkPoints> points = load_points(field_vertices);
  vtkNew<vtkPoints> linePoints;

  std::vector<std::array<double, 3>> b(field_vertices.size());
  std::vector<double> weights(field_vertices.size());
  for (int i = 0; i < field_vertices.size(); i++) {
    b[i] = field_vertices[i] + field_normals[i];
    weights[i] = (dot(field_normals[i], {0.3, 0.3, 0.3}) + 1) / 2;
  };

  vtkNew<vtkCellArray> linesArray = load_lines(field_vertices, b, linePoints);
  vtkNew<vtkDoubleArray> lineScalars = load_scalars(weights);
  lineScalars->SetNumberOfComponents(1);

  // Render
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  render_points(points, renderer);
  render_lines(linesArray, lineScalars, linePoints, renderer);

  // Window Render
  renderWindow->Render();
  renderWindowInteractor->Start();

  return 0;
}
