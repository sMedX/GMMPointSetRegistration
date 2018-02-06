#include <itkMesh.h>
#include <itkTransformMeshFilter.h>
#include <itkLBFGSOptimizer.h>

#include "itkGMMPointSetToPointSetRegistrationMethod.h"
#include "itkPointSetPropertiesCalculator.h"
#include "itkInitializeTransform.h"
#include "itkInitializeMetric.h"
#include "itkPointSetToPointSetMetrics.h"

#include "itkIOutils.h"
#include "argsCustomParsers.h"

const unsigned int Dimension = 3;
typedef itk::Mesh<float, Dimension> FixedMeshType;
typedef itk::PointSet<FixedMeshType::PixelType, Dimension> FixedPointSetType;
typedef itk::Mesh<float, Dimension> MovingMeshType;
typedef itk::PointSet<MovingMeshType::PixelType, Dimension> MovingPointSetType;
typedef itk::Transform <double, Dimension, Dimension> TransformType;

int main(int argc, char** argv) {

  // parse input arguments
  args::ArgumentParser parser("GMM-based point set to point set registration.", "");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  args::Group allRequired(parser, "Required arguments:", args::Group::Validators::All);

  args::ValueFlag<std::string> argFixedFileName(allRequired, "fixed", "The fixed mesh (point-set) filename", {'f', "fixed"});
  args::ValueFlag<std::string> argMovingFileName(allRequired, "moving", "The moving mesh (point-set) filename", {'m', "moving"});
  args::ValueFlag<std::string> argOutputFileName(parser, "output", "The output mesh (point-set) filename", {'o', "output"});

  args::ValueFlag<std::vector<double>, args::DoubleVectorReader> argScale(allRequired, "scale", "The scale levels", {"scale"});
  args::ValueFlag<size_t> argNumberOfIterations(parser, "iterations", "The number of iterations", {"iterations"}, 1000);
  args::Flag trace(parser, "trace", "Optimizer iterations tracing", {"trace"});

  const std::string transformDescription =
    "The type of transform (That is number):\n"
    "  0 : Translation\n"
    "  1 : Versor3D\n"
    "  2 : Similarity\n"
    "  3 : ScaleSkewVersor3D\n";

  args::ValueFlag<size_t> argTypeOfTransform(parser, "transform", transformDescription, {'t', "transform"}, 0);
  
  const std::string metricDescription =
    "The type of metric (That is number):\n"
    "  0 : L2Rigid\n"
    "  1 : L2\n"
    "  2 : KC\n";

  args::ValueFlag<size_t> argTypeOfMetric(parser, "metric", metricDescription, {'M', "metric"}, 0);

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

  std::string fixedFileName = args::get(argFixedFileName);
  std::string movingFileName = args::get(argMovingFileName);
  size_t numberOfIterations = args::get(argNumberOfIterations);
  size_t typeOfTransform = args::get(argTypeOfTransform);
  size_t typeOfMetric = args::get(argTypeOfMetric);

  std::cout << "options" << std::endl;
  std::cout << "number of iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // read meshes
  FixedMeshType::Pointer fixedMesh = FixedMeshType::New();
  if (!readMesh<FixedMeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << fixedFileName << std::endl;
  std::cout << "number of points " << fixedMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  MovingMeshType::Pointer movingMesh = MovingMeshType::New();
  if (!readMesh<MovingMeshType>(movingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << movingFileName << std::endl;
  std::cout << "number of points " << movingMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  FixedPointSetType::Pointer fixedPointSet = FixedPointSetType::New();
  fixedPointSet->SetPoints(fixedMesh->GetPoints());

  MovingPointSetType::Pointer movingPointSet = MovingPointSetType::New();
  movingPointSet->SetPoints(movingMesh->GetPoints());

  //--------------------------------------------------------------------
  // initialize scales
  typedef itk::PointSetPropertiesCalculator<FixedPointSetType> FixedPointSetPropertiesCalculatorType;
  FixedPointSetPropertiesCalculatorType::Pointer fixedPointSetCalculator = FixedPointSetPropertiesCalculatorType::New();
  fixedPointSetCalculator->SetPointSet(fixedPointSet);
  fixedPointSetCalculator->Compute();
  fixedPointSetCalculator->PrintReport(std::cout);

  typedef itk::PointSetPropertiesCalculator<MovingPointSetType> MovingPointSetPropertiesCalculatorType;
  MovingPointSetPropertiesCalculatorType::Pointer movingPointSetCalculator = MovingPointSetPropertiesCalculatorType::New();
  movingPointSetCalculator->SetPointSet(movingPointSet);
  movingPointSetCalculator->Compute();
  movingPointSetCalculator->PrintReport(std::cout);

  itk::Array<double> scale(args::get(argScale).size());
  for (size_t n = 0; n < scale.size(); ++n) {
    scale[n] = args::get(argScale)[n] * movingPointSetCalculator->GetScale();
  }

  // initialize transform
  typedef itk::VersorRigid3DTransform<double> InitialTransformType;
  InitialTransformType::Pointer fixedInitialTransform = InitialTransformType::New();
  fixedInitialTransform->SetCenter(fixedPointSetCalculator->GetCenter());
  fixedInitialTransform->SetIdentity();

  InitialTransformType::Pointer movingInitialTransform = InitialTransformType::New();
  movingInitialTransform->SetCenter(movingPointSetCalculator->GetCenter());
  movingInitialTransform->SetIdentity();

  typedef itk::InitializeTransform<double> TransformInitializerType;
  TransformInitializerType::Pointer transformInitializer = TransformInitializerType::New();
  transformInitializer->SetMovingLandmark(movingPointSetCalculator->GetCenter());
  transformInitializer->SetFixedLandmark(fixedPointSetCalculator->GetCenter());
  transformInitializer->SetTypeOfTransform(typeOfTransform);
  transformInitializer->Update();
  transformInitializer->PrintReport();
  TransformType::Pointer transform = transformInitializer->GetTransform();
  std::cout << " fixed " << fixedPointSetCalculator->GetCenter() << std::endl;
  std::cout << "moving " << movingPointSetCalculator->GetCenter() << std::endl;
  std::cout << " scale " << scale << std::endl;
  //--------------------------------------------------------------------
  // initialize optimizer
  typedef itk::LBFGSOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMaximumNumberOfFunctionEvaluations(numberOfIterations);
  optimizer->SetScales(transformInitializer->GetScales());
  optimizer->SetTrace(trace);
  optimizer->MinimizeOn();

  //--------------------------------------------------------------------
  // metric
  typedef itk::InitializeMetric<FixedPointSetType, MovingPointSetType> InitializeMetricType;
  InitializeMetricType::Pointer metricInitializer = InitializeMetricType::New();
  metricInitializer->SetTypeOfMetric(typeOfMetric);
  try {
    metricInitializer->Initialize();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  metricInitializer->PrintReport();
  //--------------------------------------------------------------------
  // perform registration
  typedef itk::GMMPointSetToPointSetRegistrationMethod<FixedPointSetType, MovingPointSetType> GMMPointSetToPointSetRegistrationMethodType;
  GMMPointSetToPointSetRegistrationMethodType::Pointer registration = GMMPointSetToPointSetRegistrationMethodType::New();
  registration->SetFixedPointSet(fixedPointSet);
  registration->SetFixedInitialTransform(fixedInitialTransform);
  registration->SetMovingPointSet(movingPointSet);
  registration->SetMovingInitialTransform(movingInitialTransform);
  registration->SetScale(scale);
  registration->SetOptimizer(optimizer);
  registration->SetMetric(metricInitializer->GetMetric());
  registration->SetTransform(transform);
  try {
    registration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << std::endl;
  std::cout << "optimizer " << registration->GetOptimizer()->GetNameOfClass() << std::endl;
  std::cout << "   scales " << registration->GetOptimizer()->GetScales() << std::endl;
  std::cout << std::endl;
  std::cout << "transform " << registration->GetTransform()->GetNameOfClass() << std::endl;
  std::cout << "Initial transform parameters " << registration->GetInitialTransformParameters() << std::endl;
  std::cout << "  Final transform parameters " << registration->GetFinalTransformParameters() << std::endl;
  std::cout << std::endl;
  std::cout << "metric " << registration->GetMetric()->GetNameOfClass() << std::endl;
  std::cout << "   Initial metric values " << registration->GetInitialMetricValues() << std::endl;
  std::cout << "     Final metric values " << registration->GetFinalMetricValues() << std::endl;
  std::cout << std::endl;

  // transform moving mesh
  typedef itk::TransformMeshFilter<MovingMeshType, MovingMeshType, TransformType> TransformMeshFilterType;
  TransformMeshFilterType::Pointer transformMesh = TransformMeshFilterType::New();
  transformMesh->SetInput(movingMesh);
  transformMesh->SetTransform(transform);
  try {
    transformMesh->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  if (argOutputFileName) {
    std::string fileName = args::get(argOutputFileName);
    std::cout << "write output mesh to the file " << fileName << std::endl;
    std::cout << std::endl;

    if (!writeMesh<MovingMeshType>(transformMesh->GetOutput(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  // compute metrics
  typedef itk::PointSetToPointSetMetrics<FixedPointSetType, MovingPointSetType> PointSetToPointSetMetricsType;
  PointSetToPointSetMetricsType::Pointer metrics = PointSetToPointSetMetricsType::New();
  metrics->SetFixedPointSet(fixedPointSet);
  metrics->SetMovingPointSet(movingPointSet);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  metrics->SetFixedPointSet(fixedPointSet);
  metrics->SetMovingPointSet(transformMesh->GetOutput());
  metrics->Compute();
  metrics->PrintReport(std::cout);

  return EXIT_SUCCESS;
}
