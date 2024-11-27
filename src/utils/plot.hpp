#include <array>
#include <vector>
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

vtkNew<vtkPoints> load_points(const std::vector<std::array<double, 3>> &points);

vtkNew<vtkDoubleArray> load_scalars(const std::vector<double> &weights);

vtkNew<vtkCellArray> load_lines(const std::vector<std::array<double, 3>> &a,
                                const std::vector<std::array<double, 3>> &b,
                                vtkNew<vtkPoints> &points);
