#include <itkMesh.h>
#include <itkTransformMeshFilter.h>
#include <itkCompositeTransform.h>
#include <itkEuclideanDistancePointMetric.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include <itkPointSetToPointSetRegistrationMethod.h>

#include "itkPointSetPropertiesCalculator.h"
#include "itkInitializeTransform.h"
#include "itkPointSetToPointSetMetrics.h"

#include "itkIOutils.h"
#include "argsCustomParsers.h"

const unsigned int Dimension = 3;
typedef itk::Mesh<float, Dimension> MeshType;
typedef itk::PointSet<MeshType::PixelType, MeshType::PointDimension> PointSetType;
typedef itk::Transform <double, MeshType::PointDimension, MeshType::PointDimension> TransformType;

int main(int argc, char** argv) {
  args::ArgumentParser parser("ICP PointSet Registration", "");
  args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });

  args::Group allRequired(parser, "Required arguments:", args::Group::Validators::All);

  args::ValueFlag<std::string> argFixedFileName(allRequired, "fixed", "The fixed mesh (point-set) filename", { 'f', "fixed" });
  args::ValueFlag<std::string> argMovingFileName(allRequired, "moving", "The moving mesh (point-set) filename", { 'm', "moving" });
  args::ValueFlag<std::string> argOutputFileName(parser, "output", "The output mesh (point-set) filename", { 'o', "output" });

  args::ValueFlag<size_t> argNumberOfIterations(parser, "iterations", "The number of iterations", { 'i', "iterations" }, 1000);
  args::Flag trace(parser, "trace", "Optimizer iterations tracing", {"trace"});

  const std::string transformDescription =
    "The type of transform (That is number):\n"
    "  0 : Translation\n"
    "  1 : Versor3D\n"
    "  2 : Similarity\n"
    "  3 : ScaleSkewVersor3D\n";

  args::ValueFlag<size_t> argTypeOfTransform(parser, "transform", transformDescription, { 't', "transform" }, 0);

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

  std::cout << "options" << std::endl;
  std::cout << "number of iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  //
  MeshType::Pointer fixedMesh = MeshType::New();
  if (!readMesh<MeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << fixedFileName << std::endl;
  std::cout << "number of points " << fixedMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  MeshType::Pointer movingMesh = MeshType::New();
  if (!readMesh<MeshType>(movingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }

  std::cout << movingFileName << std::endl;
  std::cout << "number of points " << movingMesh->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //--------------------------------------------------------------------
  // read meshes
  if (!readMesh<MeshType>(movingMesh, movingFileName)) {
    return EXIT_FAILURE;
  }
  PointSetType::Pointer movingPointSet = PointSetType::New();
  movingPointSet->SetPoints(movingMesh->GetPoints());

  if (!readMesh<MeshType>(fixedMesh, fixedFileName)) {
    return EXIT_FAILURE;
  }
  PointSetType::Pointer fixedPointSet = PointSetType::New();
  fixedPointSet->SetPoints(fixedMesh->GetPoints());

  // initialize transform
  typedef itk::PointSetPropertiesCalculator<PointSetType> PointSetPropertiesCalculatorType;
  PointSetPropertiesCalculatorType::Pointer fixedPointSetCalculator = PointSetPropertiesCalculatorType::New();
  fixedPointSetCalculator->SetPointSet(fixedPointSet);
  fixedPointSetCalculator->Compute();
  fixedPointSetCalculator->Print(std::cout);

  PointSetPropertiesCalculatorType::Pointer movingPointSetCalculator = PointSetPropertiesCalculatorType::New();
  movingPointSetCalculator->SetPointSet(movingPointSet);
  movingPointSetCalculator->Compute();
  movingPointSetCalculator->Print(std::cout);

  typedef itk::Similarity3DTransform<double> InitialTransformType;
  InitialTransformType::Pointer fixedInitialTransform;
  InitialTransformType::Pointer movingInitialTransform;

  // initialize transform
  typedef itk::InitializeTransform<double> TransformInitializerType;
  TransformInitializerType::Pointer transformInitializer = TransformInitializerType::New();
  transformInitializer->SetMovingLandmark(movingPointSetCalculator->GetCenter());
  transformInitializer->SetFixedLandmark(fixedPointSetCalculator->GetCenter());
  transformInitializer->SetTypeOfTransform(typeOfTransform);
  transformInitializer->Initialize();
  transformInitializer->Print(std::cout);
  TransformType::Pointer transform = transformInitializer->GetTransform();

  //--------------------------------------------------------------------
  // initialize optimizer
  typedef itk::LevenbergMarquardtOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetUseCostFunctionGradient(false);

  OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
  scales.Fill(0.01);

  const double        gradientTolerance = 1e-5;    // convergence criterion
  const double        valueTolerance = 1e-5;    // convergence criterion
  const double        epsilonFunction = 1e-6;   // convergence criterion

  optimizer->SetScales(scales);
  optimizer->SetNumberOfIterations(numberOfIterations);
  optimizer->SetValueTolerance(valueTolerance);
  optimizer->SetGradientTolerance(gradientTolerance);
  optimizer->SetEpsilonFunction(epsilonFunction);
  optimizer->SetDebug(trace);

  //--------------------------------------------------------------------
  // metric
  typedef itk::EuclideanDistancePointMetric<PointSetType, PointSetType> MetricType;
  MetricType::Pointer  metric = MetricType::New();
  metric->SetFixedPointSet(fixedPointSet);
  metric->SetMovingPointSet(movingPointSet);

  //--------------------------------------------------------------------
  // perform registration
  typedef itk::PointSetToPointSetRegistrationMethod<PointSetType, PointSetType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();
  registration->SetInitialTransformParameters(transform->GetParameters());
  registration->SetFixedPointSet(fixedPointSet);
  registration->SetMovingPointSet(movingPointSet);
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform(transform);
  try {
    registration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }
  registration->GetTransform();

  std::cout << std::endl;
  std::cout << registration->GetMetric()->GetNameOfClass() << std::endl;
  std::cout << optimizer->GetStopConditionDescription() << std::endl;
  std::cout << "initial parameters " << registration->GetInitialTransformParameters() << std::endl;
  std::cout << "  final parameters " << optimizer->GetCurrentPosition() << std::endl;
  std::cout << "             value " << optimizer->GetValue() << std::endl;
  std::cout << std::endl;

  typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformMeshFilterType;
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

    if (!writeMesh<MeshType>(transformMesh->GetOutput(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  PointSetType::Pointer outputPointSet = PointSetType::New();
  outputPointSet->SetPoints(transformMesh->GetOutput()->GetPoints());

  typedef itk::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
  PointSetToPointSetMetricsType::Pointer metrics = PointSetToPointSetMetricsType::New();
  metrics->SetFixedPointSet(fixedPointSet);
  metrics->SetMovingPointSet(movingPointSet);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  metrics->SetFixedPointSet(fixedPointSet);
  metrics->SetMovingPointSet(outputPointSet);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  return EXIT_SUCCESS;
}
