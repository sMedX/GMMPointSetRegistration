#include <itkMesh.h>
#include <itkEuler3DTransform.h>
#include <itkTransformMeshFilter.h>

#include "args.hxx"
#include "argsCustomParsers.h"
#include "itkIOutils.h"

typedef itk::Mesh<float, 3U> MeshType;
typedef itk::Euler3DTransform <double> TransformType;

int main(int argc, char** argv) {

  args::ArgumentParser parser("Application to transform mesh", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Group allRequired(parser, "Required arguments:", args::Group::Validators::All);
  args::ValueFlag<std::string> argInputFile(allRequired, "input", "The input mesh file name", {'i', "input"});
  args::ValueFlag<std::string> argOutputFile(allRequired, "output", "The output mesh file name", {'o', "output"});
  args::ValueFlag<std::vector<double>, args::DoubleVectorReader> argParameters(allRequired, "parameters", "The parameters to transform mesh", {'p', "parameters"});

  try {
    parser.ParseCLI(argc, argv);
  }
  catch (args::Help) {
    std::cout << parser;
    return EXIT_SUCCESS;
  }
  catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return EXIT_FAILURE;
  }
  catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return EXIT_FAILURE;
  }

  std::string inputFile = args::get(argInputFile);
  std::string outputFile = args::get(argOutputFile);

  MeshType::Pointer mesh = MeshType::New();
  if (!readMesh<MeshType>(mesh, inputFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input mesh " << inputFile << std::endl;
  std::cout << "number of points " << mesh->GetNumberOfPoints() << std::endl;
  std::cout << "number of cells  " << mesh->GetNumberOfCells() << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // initialize transform
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  transform->SetCenter(mesh->GetBoundingBox()->GetCenter());

  if (args::get(argParameters).size() < transform->GetNumberOfParameters()) {
    std::cout << "The number of parameters to transform input mesh must be equal to " << transform->GetNumberOfParameters() << std::endl;
    return EXIT_FAILURE;
  }

  TransformType::ParametersType parameters;
  parameters.set_size(transform->GetNumberOfParameters());

  for (int n = 0; n < transform->GetNumberOfParameters(); ++n) {
    parameters[n] = args::get(argParameters)[n];
  }

  transform->SetParameters(parameters);

  std::cout << "transform  " << transform->GetNameOfClass() << std::endl;
  std::cout << "center     " << transform->GetCenter() << std::endl;
  std::cout << "parameters " << transform->GetParameters() << std::endl;
  std::cout << std::endl;

  typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformMeshFilterType;
  TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
  transformMeshFilter->SetInput(mesh);
  transformMeshFilter->SetTransform(transform);
  try {
    transformMeshFilter->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "output mesh " << outputFile << std::endl;
  std::cout << std::endl;
  if (!writeMesh<MeshType>(transformMeshFilter->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
