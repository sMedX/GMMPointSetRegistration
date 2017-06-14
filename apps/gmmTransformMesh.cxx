#include <itkMesh.h>
#include <itkEuler3DTransform.h>
#include <itkTransformMeshFilter.h>

#include "args.hxx"
#include "agtkIO.h"

typedef itk::Mesh<float, 3U> MeshType;
typedef itk::Euler3DTransform <double> TransformType;

int main(int argc, char** argv) {

  args::ArgumentParser parser("GMM Transform Mesh", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  
  args::Group allRequired(parser, "Required arguments:", args::Group::Validators::All);

  args::ValueFlag<std::string> argInputFile(allRequired, "input", "The input mesh filename", {'i', "input"});
  args::ValueFlag<std::string> argOutputFile(allRequired, "output", "The output mesh filename", {'o', "output"});
  args::ValueFlagList<double> argTransform(allRequired, "transform", "The transform vector", {'t', "transform"});

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
  std::vector<double> options = args::get(argTransform);

  MeshType::Pointer model = MeshType::New();
  if (!agtk::readMesh<MeshType>(model, inputFile)) {
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
  if (!agtk::writeMesh<MeshType>(transformMeshFilter->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
