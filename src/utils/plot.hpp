#include <array>
#include <vector>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3DMapper.h>
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
vtkSmartPointer<vtkPoints> combine_points(vtkPoints *p1, vtkPoints *p2);

vtkNew<vtkDoubleArray> load_scalars(const std::vector<double> &weights);
vtkSmartPointer<vtkDoubleArray> combine_scalars(vtkDoubleArray *p1,
                                                vtkDoubleArray *p2);

vtkNew<vtkCellArray> load_lines(const std::vector<std::array<double, 3>> &a,
                                const std::vector<std::array<double, 3>> &b,
                                vtkNew<vtkPoints> &points);
