
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <time.h>
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
			int offset = floor(pow(normalRandom() * spread, 2) + 0.5);
			offset += nails; // to eliminate negative offsets and %
			dna[i] = (dna[i] + offset) % nails;
		}
	}

	void draw(Image &canvas, const std::vector<Point> &nails, const float threadOpacity) {
		canvas.fill(0.0f);
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
				float c = canvas(x0, y0, 0, 0);
				c += threadOpacity;
				if (c > 1.0f) c = 1.0f;
				canvas(x0, y0, 0, 0) = c;				

				if ((x0 == x1) && (y0 == y1)) break;
				int e2 = 2 * err;
				if (e2 >-dy) { err -= dy; x0 += sx; }
				if (e2 < dx) { err += dx; y0 += sy; }
			}
		}
		//canvas.normalize(0.0, 1.0);
		//canvas.display();
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

	double calculateFitness(const Image &kernel, Image &canvas, const std::vector<Point> &nails, const float threadOpacity) {
		double minFitness = 10e+7;
		double maxFitness = -10e+7;
		double average = 0.0;

		for (int i = 0; i < population.size(); i++) {
			//printf("%i\n", i);
			population[i].draw(canvas, nails, threadOpacity);
			double fitness = 0.0;
			double error = 0.0;
			for (int x = 0; x < kernel.width(); x++) {
				for (int y = 0; y < kernel.height(); y++) {
					//fitness += kernel(x, y) * canvas(x, y);
					fitness += (kernel(x, y, 0, 0) - 0.5) * (canvas(x, y, 0, 0) - 0.5);
					//error += pow(kernel(x, y, 0, 0) - canvas(x, y, 0, 0), 2);
				}
			}
			//fitness = 100.0 / (error + 1);
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

		//printf("Max: %.2lf\tMin: %.2lf\n", maxFitness, minFitness);

		return average / population.size();
	}

	Genome* getParent() {
		double random = uniformRandom() * weightSum;
		int i;
		for (i = 0; i < population.size() - 1; i++) {
			if (random > population[i].weight)
				break;
		}
		return &(population[i]);
	}
};

int getCrossoverJumpSize(int size) {
	return floor(uniformRandom() * size / 5 + size / 5);
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

void generateNewPopulation(Population &base, Population &result, int nails) {
	for (int i = 0; i < result.population.size(); i++) {
		// select 2 parents
		Genome* mother = base.getParent();
		Genome* father = base.getParent();

		// crossover
		result.population[i] = crossover(mother, father);

		// mutate
		result.population[i].mutate(2.0, nails);
	}
}

void run() {
	int nails = 150;
	int populationSize = 200;
	int genomeSize = nails * 4;
	int imageSize = 100;
	float threadOpacity = 1.0 / 10;

	char *filename = "stripes.png";

	CImg<unsigned char> file(filename);
	Image kernel(imageSize, imageSize, 1, 1, 0.0f);

	file.resize(imageSize, imageSize);

	for (int i = 0; i < imageSize; i++) {
		for (int j = 0; j < imageSize; j++) {
			double c = 0.299 * file(i, j, 0, 0) + 0.587 * file(i, j, 0, 1) + 0.114 * file(i, j, 0, 2);
			//c = 1.0 - c / 255.0;
			c = 1.0 - c / 128.0;
			kernel(i, j, 0, 0) = c;

			int dx = i - imageSize / 2;
			int dy = j - imageSize / 2;
			if (dx * dx + dy * dy >= (imageSize * imageSize / 4) * 0.8) {
				kernel(i, j, 0, 0) = 0.0;
			}
		}
	}

	Image canvas(imageSize, imageSize, 1, 1, 0.0f);

	std::vector<Point> nailPoints;
	nailPoints.resize(nails);
	for (int i = 0; i < nails; i++) {
		double a = 3.14159265359 * 2 * i / nails;
		double s = std::sin(a);
		double c = std::cos(a);
		double r = imageSize / 2;
		nailPoints[i].first = floor(r + c * r);
		nailPoints[i].second = floor(r + s * r);
	}

	Population populations[2];
	populations[0].init(populationSize, genomeSize, nails);
	populations[1].init(populationSize, genomeSize, nails);

	double fitness;// = population.calculateFitness(kernel, canvas, nailPoints);
	int generation = 0;
	
	while (generation < 2020) {
		int current = generation % 2;
		int next = (generation + 1) % 2;
		fitness = populations[current].calculateFitness(kernel, canvas, nailPoints, threadOpacity);

		printf("Generation: %i,\tFitness: %.7lf,\tMax: %.7lf\tMin: %.7lf\n", generation, fitness, populations[current].population.front().fitness, populations[current].population.back().fitness);

		generateNewPopulation(populations[current], populations[next], nails);

		generation++;
	}
	fitness = populations[generation % 2].calculateFitness(kernel, canvas, nailPoints, threadOpacity);
	printf("Generation: %i,\tFitness: %lf\n", generation, fitness);

	Genome *win = &(populations[generation % 2].population.front());
	win->draw(canvas, nailPoints, threadOpacity);
	canvas.normalize(0, 255);
	canvas *= -1;
	canvas += 255;
	canvas.display("Lowest");

	win = &(populations[generation % 2].population.back());
	win->draw(canvas, nailPoints, threadOpacity);
	canvas.normalize(0, 255);
	canvas *= -1;
	canvas += 255;
	canvas.display("Highest");
}

int main(int argc, char **argv) {
	srand(clock());
	/*for (int i = 0; i < 20; i++) {
		printf("%i\n", (int)floor(normalRandom() * 2 + 0.5));
	}*/

	run();

	char c;
	scanf_s("%c", &c);
}