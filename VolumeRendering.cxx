#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkImageData.h>
#include <vtkCellLocator.h>
#include <vtkRectilinearGridReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkProbeFilter.h>
#include <vtkBox.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPNGWriter.h>
#include <Kokkos_Core.hpp>
#include <fstream>

#include "DataStructs.h"

using namespace std;

const double PI = 3.1415926;

double colorEq(double front_o, double back_o, double front, double back) {
	return front + (1 - front_o) * back_o * back;
	// 0.0000092 + (.002355) *  0.304082 * 0
}

double map225to1(unsigned char val) {
	int v = val;
	double slope = 1.0 / 255.0;
	/*cout << "v: " << v << endl;
	cout << "double v: " << (double)v << endl;
	cout << "mapped color " << v << " to " << slope * (double) v << endl;*/
	return slope * (double) v;
}

unsigned char map1to255(double val) {
	double slope = 255.0;
	return (unsigned char)(slope * val);
}

void sampleAlongRay(Vector3<double> ray, const int n_samples, 
					Camera cam, int step_size, double* samples,
					vtkRectilinearGrid* data, vtkCellLocator* locator,
					vtkDataArray* scalars) {
	for (int i = 0; i < n_samples; i++) {
		Vector3<double> point = cam.position + (ray * (cam.near + (double)i * step_size));
		//cout << point << endl;
		vtkIdType cellId = locator->FindCell(point.coords);
		//cout << cellId << endl;

		vtkCell* cell = data->GetCell(cellId);
		double closest_pt[3];
		double pcoords[3];
		double weights[8];
		double dist2;
		int subId;
		int is_cell = cell->EvaluatePosition(point.coords, closest_pt, subId, pcoords, dist2, weights);
		double interp_val = 0.0;

		if (is_cell) {
			for (int j = 0; j < 8; j++) {
				vtkIdType p = cell->GetPointId(j);
				double val = scalars->GetTuple1(p);
				interp_val += weights[j] * val;
			}
		}
		else {
			interp_val = 0.0;
		}
		samples[i] = interp_val;
	}
}

int main(int argc, char* argv[])
{
	Kokkos::initialize(argc, argv);
	{
		// Allocate a 1-dimensional view of integers
		Kokkos::View<int*> v("v", 5);
		// Fill view with sequentially increasing values v=[0,1,2,3,4]
		Kokkos::parallel_for("fill", 5, KOKKOS_LAMBDA(int i) { v(i) = i; });
		// Compute accumulated sum of v's elements r=0+1+2+3+4
		int r;
		Kokkos::parallel_reduce(
			"accumulate", 5,
			KOKKOS_LAMBDA(int i, int& partial_r) { partial_r += v(i); }, r);
		// Check the result
		KOKKOS_ASSERT(r == 10);
	}
	Kokkos::printf("Goodbye World\n");
	Kokkos::finalize();

	const int IMG_SIZE = 100;
	const int REF_SAMPLE = 500;
	const int SAMPLES_PER_RAY = 256;
	double op_ratio = (double)REF_SAMPLE / (double)SAMPLES_PER_RAY;
	const char* FILE_NAME = "astro64.vtk";

	auto image = vtkImageData::New();
	image->SetDimensions(IMG_SIZE, IMG_SIZE, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	unsigned char* image_buffer = (unsigned char*) image->GetScalarPointer();

	ifstream inputFile(FILE_NAME);
	Camera cam = SetupCamera();
	TransferFunction tf = SetupTransferFunction();

	// vector initialization for finding rays
	Vector3<double> look = cam.focus - cam.position;
	Vector3<double> look_norm = look / look.Magnitude();

	Vector3<double> lxu = look.Cross(cam.up);
	Vector3<double> ru = lxu / lxu.Magnitude();
	Vector3<double> lxru = look.Cross(ru);
	Vector3<double> rv = lxru / lxru.Magnitude();
	double rads = cam.angle * (PI / 180.0);
	double k = (2.0 * tan(rads / 2.0)) / (double)IMG_SIZE;
	Vector3<double> dx = k * ru;
	Vector3<double> dy = k * rv;

	// getting extents to calculate step size
	vtkSmartPointer<vtkRectilinearGridReader> reader =
		vtkSmartPointer<vtkRectilinearGridReader>::New();
	reader->SetFileName(FILE_NAME); reader->Update();
	vtkRectilinearGrid* data = reader->GetOutput();
	double bounds[6];
	data->GetBounds(bounds);


	double step_size = (cam.far - cam.near) / (double)(SAMPLES_PER_RAY - 1);
	cout << "step size: " << step_size << endl;

	cout << "ru: " << ru << endl;
	cout << "rv: " << rv << endl;
	cout << "dx " << dx << endl;
	cout << "dy: " << dy << endl;

	vtkCellLocator* locator = vtkCellLocator::New();
	locator->SetDataSet(data);


	for (int y = 0; y < IMG_SIZE; y++) {
		for (int x = 0; x < IMG_SIZE; x++) {
			Vector3<double> ray = look_norm + ((2.0 * (double) x + 1.0 - (double) IMG_SIZE) / 2.0) * dx
											+ ((2.0 * (double) y + 1.0 - (double) IMG_SIZE) / 2.0) * dy;
			Vector3<double> end_pt(cam.position.x + ray.x * 1e8,
								   cam.position.y + ray.y * 1e8,
							   	   cam.position.z + ray.z * 1e8);

			double sample[SAMPLES_PER_RAY];
			vtkPointData* pt_data = data->GetPointData();
			string arr_name = pt_data->GetArrayName(0);
			vtkDataArray* scalars = pt_data->GetArray(arr_name.c_str());


			sampleAlongRay(ray, SAMPLES_PER_RAY, cam, step_size, sample, data, locator, scalars);

			// calculate color
			double sampleRGB[3]; 
			double final_opacity;

			for (int i = 0; i < SAMPLES_PER_RAY; i++) {
				unsigned char RGB[3];
				double opacity;
				tf.ApplyTransferFunction(sample[i], RGB, opacity);

				if (i >= 1) {
					opacity = 1 - pow(1 - opacity, op_ratio);
					sampleRGB[0] = colorEq(final_opacity, opacity, sampleRGB[0], map225to1(RGB[0]));
					sampleRGB[1] = colorEq(final_opacity, opacity, sampleRGB[1], map225to1(RGB[1]));
					sampleRGB[2] = colorEq(final_opacity, opacity, sampleRGB[2], map225to1(RGB[2]));
					final_opacity = final_opacity + (1 - final_opacity) * opacity;
				}
				else {
					sampleRGB[0] = map225to1(RGB[0]);
					sampleRGB[1] = map225to1(RGB[1]);
					sampleRGB[2] = map225to1(RGB[2]);
					final_opacity = opacity;
				}
			}

			int img_index = (y * IMG_SIZE + x) * 3;
			image_buffer[img_index] = map1to255(sampleRGB[0]);
			image_buffer[img_index+1] = map1to255(sampleRGB[1]);
			image_buffer[img_index+2] = map1to255(sampleRGB[2]);

		}
	}

	auto writer = vtkPNGWriter::New();
	writer->SetFileName("test_out.png");
	writer->SetInputData(image);
	writer->Write();
	
}


