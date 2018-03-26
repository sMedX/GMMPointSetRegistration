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

  args::Group required(parser, "Required arguments:", args::Group::Validators::All);

  args::ValueFlag<std::string> argFixedFileName(required, "fixed", "The fixed mesh (point-set) file name", {'f', "fixed"});
  args::ValueFlag<std::string> argMovingFileName(required, "moving", "The moving mesh (point-set) file name", {'m', "moving"});
  args::ValueFlag<std::string> argOutputFileName(parser, "output", "The output mesh (point-set) file name", {'o', "output"});
  args::ValueFlag<std::string> argTargetFileName(parser, "target", "The target mesh (point-set) file name", {'t', "target"});

  args::ValueFlag<std::vector<double>, args::DoubleVectorReader> argScale(required, "scale", "The scale levels", {"scale"});
  args::ValueFlag<size_t> argNumberOfIterations(parser, "iterations", "The number of iterations", {"iterations"}, 1000);
  args::Flag trace(parser, "trace", "Optimizer iterations tracing", {"trace"});

  const std::string transformDescription =
    "The type of transform (That is number):\n"
    "  0 : Translation\n"
    "  1 : Versor3D\n"
    "  2 : Similarity\n"
    "  3 : ScaleSkewVersor3D\n";

  args::ValueFlag<size_t> argTypeOfTransform(parser, "transform", transformDescription, {"transform"}, 0);
  
  const std::string metricDescription =
    "The type of metric (That is number):\n"
    "  0 : L2Rigid\n"
    "  1 : L2\n"
    "  2 : KC\n";

  args::ValueFlag<size_t> argTypeOfMetric(parser, "metric", metricDescription, {"metric"}, 0);

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
  std::string targetFileName = args::get(argTargetFileName);
  size_t numberOfIterations = args::get(argNumberOfIterations);
  size_t typeOfTransform = args::get(argTypeOfTransform);
  size_t typeOfMetric = args::get(argTypeOfMetric);

  std::cout << "options" << std::endl;
  std::cout << "number of iterations " << numberOfIterations << std::endl;
  std::cout << "transform " << typeOfTransform << std::endl;
  std::cout << "metric    " << typeOfMetric << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // read meshes
  FixedMeshType::Pointer fixedMesh = FixedMeshType::New();
  if (!readMesh<FixedMeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }

  FixedPointSetType::Pointer fixedPointSet = FixedPointSetType::New();
  fixedPointSet->SetPoints(fixedMesh->GetPoints());

  std::cout << "fixed mesh " << fixedFileName << std::endl;
  std::cout << "number of points " << fixedMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  MovingMeshType::Pointer movingMesh = MovingMeshType::New();
  if (!readMesh<MovingMeshType>(movingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }

  MovingPointSetType::Pointer movingPointSet = MovingPointSetType::New();
  movingPointSet->SetPoints(movingMesh->GetPoints());

  std::cout << "moving mesh " << movingFileName << std::endl;
  std::cout << "number of points " << movingMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  FixedMeshType::Pointer targetMesh = nullptr;
  FixedPointSetType::Pointer targetPointSet = nullptr;
  if (argTargetFileName) {
    targetMesh = FixedMeshType::New();
    if (!readMesh<FixedMeshType>(targetMesh, targetFileName)) {
      return EXIT_FAILURE;
    }

    targetPointSet = FixedPointSetType::New();
    targetPointSet->SetPoints(targetMesh->GetPoints());

    std::cout << "target mesh " << targetFileName << std::endl;
    std::cout << "number of points " << targetMesh->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;
  }

  //--------------------------------------------------------------------
  // initialize scales
  typedef itk::PointSetPropertiesCalculator<FixedPointSetType> FixedPointSetPropertiesCalculatorType;
  FixedPointSetPropertiesCalculatorType::Pointer fixedPointSetCalculator = FixedPointSetPropertiesCalculatorType::New();
  fixedPointSetCalculator->SetPointSet(fixedPointSet);
  try {
    fixedPointSetCalculator->Compute();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  fixedPointSetCalculator->Print(std::cout);

  typedef itk::PointSetPropertiesCalculator<MovingPointSetType> MovingPointSetPropertiesCalculatorType;
  MovingPointSetPropertiesCalculatorType::Pointer movingPointSetCalculator = MovingPointSetPropertiesCalculatorType::New();
  movingPointSetCalculator->SetPointSet(movingPointSet);
  try {
    movingPointSetCalculator->Compute();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  movingPointSetCalculator->Print(std::cout);

  itk::Array<double> scale(args::get(argScale).size());
  for (size_t n = 0; n < scale.size(); ++n) {
    scale[n] = args::get(argScale)[n] * movingPointSetCalculator->GetScale();
  }

  // initialize transform
  typedef itk::InitializeTransform<double> TransformInitializerType;
  TransformInitializerType::Pointer initializerTransform = TransformInitializerType::New();
  initializerTransform->SetMovingLandmark(movingPointSetCalculator->GetCenter());
  initializerTransform->SetFixedLandmark(fixedPointSetCalculator->GetCenter());
  initializerTransform->SetTypeOfTransform(typeOfTransform);
  try {
    initializerTransform->Initialize();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  initializerTransform->Print(std::cout);
  TransformType::Pointer transform = initializerTransform->GetTransform();

  //--------------------------------------------------------------------
  // initialize optimizer
  typedef itk::LBFGSOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMaximumNumberOfFunctionEvaluations(numberOfIterations);
  optimizer->SetScales(initializerTransform->GetScales());
  optimizer->SetTrace(trace);
  optimizer->MinimizeOn();

  //--------------------------------------------------------------------
  // metric
  typedef itk::InitializeMetric<FixedPointSetType, MovingPointSetType> InitializeMetricType;
  InitializeMetricType::Pointer initializerMetric = InitializeMetricType::New();
  initializerMetric->SetTypeOfMetric(typeOfMetric);
  try {
    initializerMetric->Initialize();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  initializerMetric->Print(std::cout);

  //--------------------------------------------------------------------
  // perform registration
  typedef itk::GMMPointSetToPointSetRegistrationMethod<FixedPointSetType, MovingPointSetType> GMMPointSetToPointSetRegistrationMethodType;
  GMMPointSetToPointSetRegistrationMethodType::Pointer registration = GMMPointSetToPointSetRegistrationMethodType::New();
  registration->SetFixedPointSet(fixedPointSet);
  registration->SetMovingPointSet(movingPointSet);
  registration->SetScale(scale);
  registration->SetOptimizer(optimizer);
  registration->SetMetric(initializerMetric->GetMetric());
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
  registration->GetMetric()->Print(std::cout);
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
  metrics->SetTargetPointSet(targetMesh);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  return EXIT_SUCCESS;
}
