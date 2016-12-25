
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string.h>
#include "CImg.h"

using namespace cimg_library;

typedef int Chromosome;
typedef std::vector<Chromosome> DNA;
typedef std::pair<int, int> Point;
typedef CImg<float> Image;

double uniformRandom()
{
	return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
}

int intRandom(int start, int end) {
	return start + rand() % (end - start);
}
// return a normally distributed random number
double normalRandom()
{
	double u1 = uniformRandom();
	double u2 = uniformRandom();
	return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

struct Triple {
	double min;
	double max;
	double mean;
};

struct Slab {
	double *data;

	Slab() {
		data = new double[width * width];
		fill(0.0);
	}

	~Slab() {
		delete[] data;
	}

	void fill(double v = 0.0) {
		for (int i = 0; i < width * width; i++) {
			data[i] = v;
		}
	}

	void fromImage(const Image &image) {
		for (int i = 0; i < image.width(); i++) {
			for (int j = 0; j < image.height(); j++) {
				if (i >= bounds[j].first && i <= bounds[j].second) {
					data[i + j * width] = image(i, j, 0, 0);
				}
			}
		}
	}

	void toImage(Image &image) {
		image.fill(0.0);
		for (int i = 0; i < width; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				image(j, i, 0, 0) = data[j + i * width];
			}
		}
	}

	void normalize() {
		double mean = 0.0;
		double dev = 0.0;
		for (int i = 0; i < width; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				mean += data[j + i * width];
			}
		}
		mean /= size;

		double d;
		for (int i = 0; i < width; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				d = data[j + i * width] - mean;
				dev += d * d;
			}
		}
		dev = sqrt(dev);

		for (int i = 0; i < width; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				data[j + i * width] = (data[j + i * width] - mean) / dev;
			}
		}
	}

	double covariate(const Slab &kernel) {
		double result = 0.0;
		for (int i = 0; i < width; i++) {
			for (int j = bounds[i].first; j <= bounds[j].second; j++) {
				result += data[j + i * width] * kernel.data[j + i * width];
			}
		}
		return result;
	}

	void addPixel(int x, int y, double opacity) {
		data[y * width + x] += opacity;
	}

	void clamp() {
		for (int i = 0; i < width; i++) {
			for (int j = bounds[i].first; j <= bounds[i].second; j++) {
				data[j + i * width] = std::min(data[j + i * width], 1.0);
			}
		}
	}

	void scan() {
		int j = width / 2;
		for (int i = 0; i < width; i++) {
			printf("%i:\t%.7lf\n", i, data[i * width + j]);
		}
	}
	
	static int width;
	static int size;
	static std::vector<Point> bounds;
	static void prepare(int imageSize) {
		width = imageSize;
		bounds.resize(width);
		int r = imageSize / 2;
		for (int i = 0; i < width; i++) {
			int start = 0;
			int end = 0;
			int flag = false;
			for (int j = 0; j < width; j++) {
				int dx = i - r;
				int dy = j - r;
				if ((dx * 2 + 1) * (dx * 2 + 1) + (dy * 2 + 1) * (dy * 2 + 1) <= (r * 2 + 1) * (r * 2 + 1)) {
					if (!flag) {
						flag = true;
						bounds[i].first = j;
					}
				} else {
					if (flag) {
						flag = false;
						bounds[i].second = j - 1;
					}
				}
			}
			if (bounds[i].second == 0)
				bounds[i].second = width - 1;
		}
		int l = 0;
		for (int i = 0; i < width; i++) {
			l += bounds[i].second - bounds[i].first;
		}
		size = l;
	}
};

int Slab::width;
int Slab::size;
std::vector<Point> Slab::bounds;

struct Genome {
	DNA dna;
	double fitness;
	double weight;

	void init(int size, int nails) {
		dna.resize(size);
		int minOffset = nails * 1 / 5;
		int maxOffset = nails * 4 / 5;
		dna[0] = rand() % nails;
		for (int i = 1; i < size; i++) {
			dna[i] = (dna[i - 1] + intRandom(minOffset, maxOffset)) % nails;
		}
		fitness = 0.0;
	}

	void mutate(double spread, int nails) {
		for (int i = 0; i < dna.size(); i++) {
			int offset = floor(pow(normalRandom() * spread * nails, 2) + 0.5);
			offset += nails; // to eliminate negative offsets and %
			dna[i] = (dna[i] + offset) % nails;
		}
	}

	void burst(float probability, int nails) {
		for (int i = 0; i < dna.size(); i++) {
			if (uniformRandom() < probability) {
				int offset = floor(uniformRandom() * nails);
				dna[i] = offset;
			}
		}
	}

	void draw(Slab &canvas, const std::vector<Point> &nails, const float threadOpacity, const std::vector<Point> &bounds) {
		canvas.fill(0.0);
		double paint = threadOpacity;
		for (int i = 0; i < dna.size() - 1; i++) {
			int x0 = nails[dna[i]].first;
			int x1 = nails[dna[i + 1]].first;
			int y0 = nails[dna[i]].second;
			int y1 = nails[dna[i + 1]].second;
			int dx = cimg::abs(x1 - x0);
			int dy = cimg::abs(y1 - y0);
			int sx = (x0 < x1) ? 1 : -1;
			int sy = (y0 < y1) ? 1 : -1;
			int err = dx - dy;
			while (true) {
				canvas.addPixel(x0, y0, paint);			

				if ((x0 == x1) && (y0 == y1)) break;
				int e2 = 2 * err;
				if (e2 >-dy) { err -= dy; x0 += sx; }
				if (e2 < dx) { err += dx; y0 += sy; }
			}
		}
		canvas.clamp();
		canvas.normalize();
	}

	// ORDER BY fitness DESC
	bool operator <(const Genome &g) {
		return fitness > g.fitness;
	}
};

typedef std::vector<Genome> Genomes;

struct Population {
	Genomes population;
	double weightSum;

	void init(int populationSize, int genomeSize, int nails) {
		population.resize(populationSize);
		for (int i = 0; i < population.size(); i++) {
			population[i].init(genomeSize, nails);
		}
	}

	double calculateFitness(const Slab &kernel, Slab &canvas, const std::vector<Point> &nails, const float threadOpacity, const std::vector<Point> &bounds) {
		double minFitness = 10e+7;
		double maxFitness = -10e+7;
		double average = 0.0;

		for (int i = 0; i < population.size(); i++) {
			//printf("%i\n", i);
			population[i].draw(canvas, nails, threadOpacity, bounds);
			double fitness = canvas.covariate(kernel);
			population[i].fitness = fitness;
			if (fitness > maxFitness) maxFitness = fitness;
			if (fitness < minFitness) minFitness = fitness;
			average += fitness;
		}
		std::sort(population.begin(), population.end());

		weightSum = 0.0;
		for (int i = 0; i < population.size(); i++) {
			double weight = (population[i].fitness - minFitness) / (maxFitness - minFitness);
			population[i].weight = weightSum;
			weightSum += weight;
		}
		
		return average / population.size();
	}

	Genome* getParent() {
		double random = uniformRandom() * weightSum / 2; // only first (successful) part, because fuck you
		int i;
		for (i = population.size() - 1; i > 0; i--) {
			if (random > population[i].weight)
				break;
		}
		//printf("Parent: %i (%lf)\n", i, population[i].fitness);
		return &(population[i]);
	}

	Genome* getParentOrdinal() {
		double random = pow(uniformRandom(), 2) * population.size() / 2;
		int i = floor(random);
		//printf("%i ", i);
		return &(population[i]);
	}

	void burst(float probability, int nails) {
		for (int i = 0; i < population.size(); i++) {
			population[i].burst(probability, nails);
		}
	}
};

double MAX_CHUNK, MIN_CHUNK;
int getCrossoverJumpSize(int size) {
	return floor(size * ( uniformRandom() * (MAX_CHUNK - MIN_CHUNK) + MIN_CHUNK) );
}

Genome crossover(Genome *mother, Genome *father) {
	Genome child;
	int size = mother->dna.size();

	child.dna.resize(size);

	const Genome *donor;
	int i = 0;
	int j = 0;
	int jump = getCrossoverJumpSize(size);
	donor = uniformRandom() < 0.5 ? mother : father;
	while (i < size) {
		child.dna[i] = donor->dna[i];
		i++;
		j++;
		if (j == jump) {
			j = 0;
			jump = getCrossoverJumpSize(size);
			donor = donor == father ? mother : father;
		}
	}
	return child;
}

void generateNewPopulation(Population &base, Population &result, int nails, int elite, float mutationSpread) {
	for (int i = 0; i < elite && i < result.population.size(); i++) {
		result.population[i] = base.population[i];
	}
	for (int i = elite; i < result.population.size(); i++) {
		// select 2 parents
		Genome* mother = base.getParentOrdinal();
		Genome* father = base.getParentOrdinal();

		// crossover
		result.population[i] = crossover(mother, father);

		// mutate
		result.population[i].mutate(mutationSpread, nails);
	}
}

Image loadImage(char *filename, int imageSize) {
	Image result(imageSize, imageSize, 1, 1, 0.0);
	CImg<unsigned char> image(filename);

	image.resize(imageSize, imageSize);

	for (int i = 0; i < imageSize; i++) {
		for (int j = 0; j < imageSize; j++) {
			double c = image(i, j, 0, 0);
			c = 1.0 - c / 255.0;
			result(i, j, 0, 0) = c;
		}
	}

	return result;
}

void run() {
	int nails = 150;
	int populationSize = 100;
	int elite = 0;
	int genomeSize = nails * 4;
	int imageSize = 100;
	double minChunk = 0.5;
	double maxChunk = 0.5;
	int iterations = 50;
	int stableIterations = 5;
	int iterationsToBurst = 100;
	float threadOpacity = 1.0 / 10;
	double mutationSpread = 1.0;
	double burstProbability = 0.005;
	double threshold = 1e-7;
	double ringDiameter = 0;
	double threadDiameter = 0;
	int saveEvery = -1;
	int historyStep = 1;
	char filename[256] = "test-i.png";

	FILE *config = NULL;
	fopen_s(&config, "config.txt", "r");
	
	char line[256];
	while (fscanf_s(config, "%s", line, 256) != EOF) {
		if (strcmp(line, "file") == 0) {
			//strcpy_s(filename, line);
			fscanf_s(config, "%s", filename, 256);
		} else
		if (strcmp(line, "nails") == 0) {
			fscanf_s(config, "%i", &nails);
		} else
		if (strcmp(line, "populationsize") == 0) {
			fscanf_s(config, "%i", &populationSize);
		} else
		if (strcmp(line, "elite") == 0) {
			fscanf_s(config, "%i", &elite);
		} else
		if (strcmp(line, "genomesize") == 0) {
			fscanf_s(config, "%i", &genomeSize);
		} else
		if (strcmp(line, "minchunk") == 0) {
			fscanf_s(config, "%lf", &minChunk);
			MIN_CHUNK = minChunk;
		} else
		if (strcmp(line, "maxChunk") == 0) {
			fscanf_s(config, "%lf", &maxChunk);
			MAX_CHUNK = maxChunk;
		} else
		if (strcmp(line, "imagesize") == 0) {
			fscanf_s(config, "%i", &imageSize);
		} else
		if (strcmp(line, "iterations") == 0) {
			fscanf_s(config, "%i", &iterations);
		}
		if (strcmp(line, "mutation") == 0) {
			fscanf_s(config, "%lf", &mutationSpread);
		} else
		if (strcmp(line, "stableiterations") == 0) {
			fscanf_s(config, "%i", &stableIterations);
		} else
		if (strcmp(line, "goodgenerations") == 0) {
			fscanf_s(config, "%i", &iterationsToBurst);
		} else
		if (strcmp(line, "threshold") == 0) {
			fscanf_s(config, "%lf", &threshold);
		} else
		if (strcmp(line, "opacity") == 0) {
			fscanf_s(config, "%lf", &threadOpacity);
		} else
		if (strcmp(line, "burst") == 0) {
			fscanf_s(config, "%lf", &burstProbability);
		}
		if (strcmp(line, "ring_diameter") == 0) {
			fscanf_s(config, "%lf", &ringDiameter);
		} else
		if (strcmp(line, "save_every") == 0) {
			fscanf_s(config, "%i", &saveEvery);
		} else
		if (strcmp(line, "history_step") == 0) {
			fscanf_s(config, "%i", &historyStep);
		} else
		if (strcmp(line, "thread_diameter") == 0) {
			fscanf_s(config, "%lf", &threadDiameter);
		} else {
			
		}
	}

	fclose(config);

	if (abs(ringDiameter) > 1e-7 && abs(threadDiameter) > 1e-7) {
		threadOpacity = imageSize * threadDiameter / ringDiameter;
	}

	std::vector<Point> nailPoints;
	nailPoints.resize(nails);
	for (int i = 0; i < nails; i++) {
		double a = 3.14159265359 * 2 * i / nails;
		double s = std::sin(a);
		double c = std::cos(a);
		double r = imageSize / 2 - 1;
		nailPoints[i].first = floor(imageSize / 2 + c * r);
		nailPoints[i].second = floor(imageSize / 2 + s * r);
	}
	std::vector<Point> bounds;
	bounds.resize(imageSize);
	for (int i = 0; i < imageSize; i++) {
		double hw = imageSize / 2 - 0.5;
		int off = sqrt(hw * hw - pow(i - hw, 2));
		bounds[i].first = floor(imageSize / 2 - off);
		bounds[i].second = floor(imageSize / 2 + off);
	}

	Population populations[2];
	populations[0].init(populationSize, genomeSize, nails);
	populations[1].init(populationSize, genomeSize, nails);

	Image kernelImage = loadImage(filename, imageSize);
	Image output(imageSize, imageSize, 1, 1, 0.0f);

	//normalizeImage(kernel);
	Slab::prepare(imageSize);

	Slab kernel;
	kernel.fromImage(kernelImage);
	kernel.normalize();

	Slab canvas;

	double fitness = 0.0;
	double prevFitness = 1.0;
	double best = -1.0;
	double prevBest = -1.0;
	int stableGeneration = 0;
	int goodGeneration = 0;
	int generation = 0;
	bool enough = false;
	std::vector<Triple> history;
	history.reserve(iterations);
	
	while (generation < iterations) {
		int current = generation % 2;
		int next = (generation + 1) % 2;

		prevFitness = fitness;
		fitness = populations[current].calculateFitness(kernel, canvas, nailPoints, threadOpacity, bounds);

		prevBest = std::max(prevBest, best);
		best = populations[current].population.front().fitness;

		Triple h;
		h.max = populations[current].population.front().fitness;
		h.min = populations[current].population.back().fitness;
		h.mean = fitness;
		history.push_back(h);
		printf("Generation: %i,\tFitness: %.7lf,\tMax: %.7lf\tMin: %.7lf\n", generation, fitness, populations[current].population.front().fitness, populations[current].population.back().fitness);
		if (saveEvery > 0 && generation % saveEvery == 0) {
			printf("SAVING\n");
			Genome *win = &(populations[generation % 2].population.front());
			win->draw(canvas, nailPoints, threadOpacity, bounds);
			canvas.toImage(output);
			char fname[100];
			sprintf_s(fname, 100, "generation-%i.bmp", generation);
			output.normalize(0, 255);
			output *= -1;
			output += 255;
			output.save(fname);
		}


		generateNewPopulation(populations[current], populations[next], nails, elite, mutationSpread);

		if (best < prevBest + 1e-10) {
			goodGeneration++;
		} else {
			goodGeneration = 0;
		}

		if (goodGeneration > iterationsToBurst) {
			goodGeneration = 0;
			populations[next].burst(burstProbability, nails);
			printf("BUSRT!\n");
		}

		if (abs(1.0 - fitness / prevFitness) < threshold) {
			stableGeneration++;
		} else {
			stableGeneration = 0;
		}
		if (stableGeneration >= stableIterations) {
			enough = true;
		}
		generation++;
	}
	fitness = populations[generation % 2].calculateFitness(kernel, canvas, nailPoints, threadOpacity, bounds);
	printf("Generation: %i,\tFitness: %lf\n", generation, fitness);

	Genome *win = &(populations[generation % 2].population.front());
	win->draw(canvas, nailPoints, threadOpacity, bounds);

	canvas.toImage(output);

	output.normalize(0, 255);
	output *= -1;
	output += 255;

	double length = 0.0;
	for (int i = 0; i < win->dna.size() - 1; i++) {
		Point p1 = nailPoints[win->dna[i]];
		Point p2 = nailPoints[win->dna[i+1]];
		int dx = p1.first - p2.first;
		int dy = p1.second - p2.second;
		length += sqrt(dx * dx + dy * dy);
	}
	length *= ringDiameter / imageSize;
	printf("Length: %.2lf m\n", length / 1000.0);

	output.save("result.bmp");

	output.display("Winner");

	FILE *outfile;
	fopen_s(&outfile, "history.csv", "w");
	for (int i = 0; i < history.size(); i++) {
		if(i % historyStep == 0)
			fprintf_s(outfile, "%i\t%.7lf\t%.7lf\t%.7lf\n", i, history[i].min, history[i].mean, history[i].max);
	}
	fclose(outfile);

	//drawHistory(history);
}

int main(int argc, char **argv) {
	srand(clock());
	/*for (int i = 0; i < 20; i++) {
		printf("%i\n", (int)floor(normalRandom() * 2 + 0.5));
	}*/

	/*Image test1 = loadImage("stripes.png", 100);
	Image test2(100, 100, 1, 1, 0.0);
	for (int i = 0; i < 100 - 3; i++) {
		for (int j = 0; j < 100; j++) {
			test2(i, j, 0, 0) = test1(i + 3, j, 0, 0);
		}
	}
	normalizeImage(test1);
	normalizeImage(test2);
	test1.display();
	test2.display();
	double r = convolve(test1, test2);
	printf("%lf", r);
	test1 -= test2;
	test1.display();*/

	run();

	char c;
	//scanf_s("%c", &c);
}