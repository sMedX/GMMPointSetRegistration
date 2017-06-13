#include <itkMesh.h>
#include <itkEuler3DTransform.h>
#include <itkTransformMeshFilter.h>

#include "utils/agtkIO.h"
#include "utils/agtkCommandLineArgumentParser.h"

using namespace agtk;

typedef itk::Mesh<float, 3U> MeshType;
typedef itk::Euler3DTransform <double> TransformType;

int main(int argc, char** argv) {

  CommandLineArgumentParser::Pointer parser = CommandLineArgumentParser::New();
  parser->SetCommandLineArguments(argc, argv);

  std::string inputFile;
  parser->GetValue("-input", inputFile);

  std::string outputFile;
  parser->GetValue("-output", outputFile);

  std::vector<double> options;
  parser->GetValue("-transform", options);

  MeshType::Pointer model = MeshType::New();
  if (!readMesh<MeshType>(model, inputFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input mesh " << inputFile << std::endl;
  std::cout << "number of points " << model->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // initialize transform
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  transform->SetCenter(model->GetBoundingBox()->GetCenter());

  TransformType::ParametersType parameters;
  for (int n = 0; n < transform->GetNumberOfParameters(); ++n) {
    parameters[n] = options[n];
  }

  transform->SetParameters(parameters);

  typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformMeshFilterType;
  TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();

  transformMeshFilter->SetInput(model);
  transformMeshFilter->SetTransform(transform);
  transformMeshFilter->Update();

  std::cout << "output mesh " << inputFile << std::endl;
  if (!writeMesh<MeshType>(transformMeshFilter->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
