#include <utility>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <map>

using namespace std;

typedef pair<string, int> quest;
typedef int id;

const int TURN_TYPE_PUSH = 0;
const int ACTOR_ID = 0;
const int VILLAIN_ID = 1;

enum class Direction {
    UP = 0,
    RIGHT = 1,
    DOWN = 2,
    LEFT = 3
};

class Dir {
public:
    static Direction Directions[4];
    static string DirectionStrings[4];

    static string GetAsString(const Direction &direction) {
        return DirectionStrings[(int) direction];
    }
};

Direction Dir::Directions[4] = {Direction::UP, Direction::RIGHT, Direction::DOWN, Direction::LEFT};
string Dir::DirectionStrings[4] = {"UP", "RIGHT", "DOWN", "LEFT"};

struct Point {
    int x{0};
    int y{0};

    Point() = default;

    Point(int x, int y) : x(x), y(y) {}

    string toString() const {
        return "{ " + to_string(x) + ", " + to_string(y) + "}";
    }

    bool operator==(const Point &other) {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Point &other) {
        return !(x == other.x && y == other.y);
    }

    static Point GetPointInDirectionOf(const Point &other, const Direction &direction) {
        switch (direction) {
            case Direction::UP:
                return PointToTopOf(other);
            case Direction::RIGHT:
                return PointToRightOf(other);
            case Direction::DOWN:
                return PointToBottomOf(other);
            case Direction::LEFT:
                return PointToLeftOf(other);
        }
    }

    static Point PointToRightOf(const Point &other) {
        return {other.x + 1, other.y};
    }

    static Point PointToLeftOf(const Point &other) {
        return {other.x - 1, other.y};
    }

    static Point PointToTopOf(const Point &other) {
        return {other.x, other.y - 1};
    }

    static Point PointToBottomOf(const Point &other) {
        return {other.x, other.y + 1};
    }
};

class Distance {
public:
    static int GetDistance(const Point &a, const Point &b) {
        return abs(a.x - b.x) + abs(a.y - b.y);
    }
};

class Tile {
public:
    static string TILE_ERR;

    static bool IsConnectedToRight(const string &other) {
        if (other == Tile::TILE_ERR) {
            return false;
        }
        bool isConnected = other[1] == '1';
        return isConnected;
    }

    static bool IsConnectedToLeft(const string &other) {
        if (other == Tile::TILE_ERR) {
            return false;
        }
        bool isConnected = other[3] == '1';
        return isConnected;
    }

    static bool IsConnectedToTop(const string &other) {
        if (other == Tile::TILE_ERR) {
            return false;
        }
        bool isConnected = other[0] == '1';
        return isConnected;
    }

    static bool IsConnectedToBottom(const string &other) {
        if (other == Tile::TILE_ERR) {
            return false;
        }
        bool isConnected = other[2] == '1';
        return isConnected;
    }
};

string Tile::TILE_ERR = "0000";

class Map {
public:
    static string GetTileAt(const string &mapStr, const Point &point) {
        int index = GetMapIndexFor(point);

        if (index == -1) {
            return Tile::TILE_ERR;
        }

        string tile = mapStr.substr(index, 4);
        return tile;
    }

    static int GetMapIndexFor(const Point &point) {
        if (point.x < 0 || point.y < 0 || point.x > 6 || point.y > 6) {
            return -1;
        }

        int index = 4 * point.x + 28 * point.y;
        return index;
    }

    static string
    PushMapInDirection(const string &map, const Direction &direction, const int id, const string &replacement,
                       string &outTile) {
        switch (direction) {
            case Direction::UP:
                return PushUpOnMap(map, id, replacement, outTile);
            case Direction::RIGHT:
                return PushRightOnMap(map, id, replacement, outTile);
            case Direction::DOWN:
                return PushDownOnMap(map, id, replacement, outTile);
            case Direction::LEFT:
                return PushLeftOnMap(map, id, replacement, outTile);
        }
    }

    static string PushRightOnMap(string map, const int id, const string &replacement, string &outTile) {
        for (int i = 6; i >= 0; --i) {
            ReplaceTileAtIndex(i, map, id, replacement, outTile, Direction::LEFT);
        }

        return map;
    }

    static string PushLeftOnMap(string map, const int id, const string &replacement, string &outTile) {
        for (int i = 0; i < 7; ++i) {
            ReplaceTileAtIndex(i, map, id, replacement, outTile, Direction::RIGHT);
        }

        return map;
    }

    static string PushDownOnMap(string map, const int id, const string &replacement, string &outTile) {
        for (int i = 6; i >= 0; --i) {
            ReplaceTileAtIndex(i, map, id, replacement, outTile, Direction::UP);
        }

        return map;
    }

    static string PushUpOnMap(string map, const int id, const string &replacement, string &outTile) {
        for (int i = 0; i < 7; ++i) {
            ReplaceTileAtIndex(i, map, id, replacement, outTile, Direction::DOWN);
        }

        return map;
    }

private:
    static void ReplaceTileAtIndex(const int index, string &map, const int id, const string &replacement,
                                   string &outTile,
                                   const Direction &replacementDirection) {

        bool isBackWard = replacementDirection == Direction::LEFT || replacementDirection == Direction::UP;
        bool isForward = !isBackWard;
        bool isVertical = replacementDirection == Direction::UP || replacementDirection == Direction::DOWN;
        bool isHorizontal = !isVertical;

        Point point = isHorizontal ?
                      /*y-axis is fixed, x moves*/ Point(index, id) :
                      /*x-axis is fixed, y moves*/ Point(id, index);


        bool isLastTile = (index == 6 && isBackWard) || (index == 0 && isForward);
        bool isFirstTile = (index == 6 && isForward) || (index == 0 && isBackWard);

        if (isLastTile) {
            outTile = GetTileAt(map, point);
        }

        int tileIndex = GetMapIndexFor(point);

        if (isFirstTile) {
            map.replace(tileIndex, 4, replacement);
        } else {
            string previousTile = GetTileAt(map, Point::GetPointInDirectionOf(point, replacementDirection));
            map.replace(tileIndex, 4, previousTile);
        }
    }
};

struct MapTile {
    const string &tile;
    const Point &point;

    MapTile(const string &tile, const Point &point) : tile(tile), point(point) {}

    bool IsConnectedToTile(const MapTile &other, Direction &outDirection) {
        if (Tile::IsConnectedToBottom(tile) && Tile::IsConnectedToTop(other.tile) &&
            Point::PointToBottomOf(point) == other.point) {
            outDirection = Direction::DOWN;
            return true;
        } else if (Tile::IsConnectedToRight(tile) && Tile::IsConnectedToLeft(other.tile) &&
                   Point::PointToRightOf(point) == other.point) {
            outDirection = Direction::RIGHT;
            return true;
        } else if (Tile::IsConnectedToLeft(tile) && Tile::IsConnectedToRight(other.tile) &&
                   Point::PointToLeftOf(point) == other.point) {
            outDirection = Direction::LEFT;
            return true;
        } else if (Tile::IsConnectedToTop(tile) && Tile::IsConnectedToBottom(other.tile) &&
                   Point::PointToTopOf(point) == other.point) {
            outDirection = Direction::UP;
            return true;
        }

        return false;
    }
};

struct Item {
    string itemName;
    Point position;
    int playerId{0};

    Item() = default;

    Item(string itemName, int x, int y, int playerId) : itemName(std::move(itemName)), position(Point(x, y)),
                                                        playerId(playerId) {}
};

struct GameState {
    string map;
    std::map<id, Point> playerPositions;
    std::map<id, string> playerTiles;
    std::map<quest, Item> items;
    std::map<id, vector<quest>> quests;

    GameState Push(const int id, const Direction &direction, const int playerId) const {
        GameState newState;
        newState.playerTiles = playerTiles;
        newState.items = items;
        newState.quests = quests;
        newState.playerPositions = playerPositions;

        //move tiles
        newState.map = Map::PushMapInDirection(map, direction, id, playerTiles.at(playerId),
                                               newState.playerTiles.at(playerId));

        bool isVerticalPush = direction == Direction::UP || direction == Direction::DOWN;

        for (const pair<quest, Item> itemPair : newState.items) {

            //Move items on pushed columns
            if (isVerticalPush && newState.items[itemPair.first].position.x == id) {
                newState.items[itemPair.first].position.y += direction == Direction::DOWN ? 1 : -1;

                //remove items that fall off map
                if (newState.items[itemPair.first].position.y > 6 ||
                    newState.items[itemPair.first].position.y < 0) {
                    newState.items[itemPair.first].position =
                            playerId == ACTOR_ID ? Point(-1, -1) : Point(-2, -2);
                }
            }
                //move items on pushed rows
            else if (!isVerticalPush && newState.items[itemPair.first].position.y == id) {
                newState.items[itemPair.first].position.x += direction == Direction::RIGHT ? 1 : -1;

                //remove items that fall off map
                if (newState.items[itemPair.first].position.x > 6 ||
                    newState.items[itemPair.first].position.x < 0) {
                    newState.items[itemPair.first].position =
                            playerId == ACTOR_ID ? Point(-1, -1) : Point(-2, -2);
                }
            }
            //move items on player replacement tiles into the map
            if (playerId == ACTOR_ID &&
                newState.items[itemPair.first].position == Point(-1, -1)) {
                newState.items[itemPair.first].position = isVerticalPush ?
                                                          Point(id, direction == Direction::DOWN ? 0 : 6) :
                                                          Point(direction == Direction::RIGHT ? 0 : 6, 1);
            }
                //move items on opponent replacement tiles into the map
            else if (playerId == VILLAIN_ID &&
                     newState.items[itemPair.first].position == Point(-2, -2)) {
                cerr << "Bringing villain item into map" << endl;
                //Bring item in
                newState.items[itemPair.first].position = isVerticalPush ?
                                                          Point(id, direction == Direction::DOWN ? 0 : 6) :
                                                          Point(direction == Direction::RIGHT ? 0 : 6, 1);
            }
        }

        for (int p_Id = 0; p_Id < 2; ++p_Id) {
            //move players on moved columns
            if (isVerticalPush && newState.playerPositions[p_Id].x == id) {
                newState.playerPositions[p_Id].y += direction == Direction::DOWN ? 1 : -1;

                //move player to other end if falls off map
                if (newState.playerPositions[p_Id].y > 6) {
                    newState.playerPositions[p_Id].y = 0;
                } else if (newState.playerPositions[p_Id].y < 0) {
                    newState.playerPositions[p_Id].y = 6;
                }
            }
                //move players on moved rows
            else if (!isVerticalPush && newState.playerPositions[p_Id].y == id) {
                newState.playerPositions[p_Id].x += direction == Direction::RIGHT ? 1 : -1;

                //move player to other end if falls off map
                if (newState.playerPositions[p_Id].x > 6) {
                    newState.playerPositions[p_Id].x = 0;
                } else if (newState.playerPositions[p_Id].y < 0) {
                    newState.playerPositions[p_Id].y = 6;
                }
            }
        }

        return newState;
    }
};

//======================================================================================================================
//  A-STAR IMPLEMENTATION
//======================================================================================================================
struct AStarNode {
    AStarNode *parent{nullptr};
    Point point;
    int distanceToGoal{0};
    int pathDistance{0};
    Direction direction{Direction::DOWN};

    AStarNode() = default;

    bool operator<(const AStarNode &rhs) {
        return totalDistance() < rhs.totalDistance();
    }

    int totalDistance() const {
        return distanceToGoal + pathDistance;
    }
};

struct AStarResult {
    vector<Direction> directions;
    bool reachedGoal{false};
    bool movedCloserToGoal{false};
    int finalDistanceToGoal{49};
    Point endPoint;
};

class AStar {
public:
    AStarResult DoSearch(const string &map, const Point &start, const Point &goal) {
        AStarResult output;
        this->DoSearch(map, start, goal, output.directions, output.reachedGoal, output.movedCloserToGoal,
                       output.finalDistanceToGoal, output.endPoint);
        return output;
    }

private:
    void DoSearch(const string &map, const Point &start, const Point &goal, vector<Direction> &outDirections,
                  bool &reachedGoal, bool &movedCloserToGoal, int &lastDistanceToGoal, Point &endPoint) const {

        std::map<string, bool> visited;

        //use greater<> to prioritize lowest
        std::priority_queue<AStarNode *, std::vector<AStarNode *>, greater<>> priorityQueue;

        reachedGoal = false;
        movedCloserToGoal = false;
        outDirections = vector<Direction>();

        AStarNode *currentNode = new AStarNode();
        currentNode->point = start;
        currentNode->distanceToGoal = Distance::GetDistance(currentNode->point, goal);

        int initialDistanceToGoal = currentNode->distanceToGoal;
        priorityQueue.push(currentNode);

        while (!priorityQueue.empty()) {
            currentNode = priorityQueue.top();
            priorityQueue.pop();

            visited[currentNode->point.toString()] = true;

            if (currentNode->point == goal) {
                reachedGoal = true;
                break;
            }

            for (auto direction: Dir::Directions) {
                string currentTile = Map::GetTileAt(map, currentNode->point);
                Point nextPoint = Point::GetPointInDirectionOf(currentNode->point, direction);
                string nextTile = Map::GetTileAt(map, nextPoint);

                if (visited.find(nextPoint.toString()) != visited.end()) {
                    continue;
                }

                Direction movementDirection = direction;
                if (MapTile(currentTile, currentNode->point)
                        .IsConnectedToTile
                                (MapTile(nextTile, nextPoint), movementDirection)) {
                    AStarNode *newNode = new AStarNode();
                    newNode->point = nextPoint;
                    newNode->parent = currentNode;
                    newNode->direction = movementDirection;
                    newNode->distanceToGoal = Distance::GetDistance(nextPoint, goal);
                    newNode->pathDistance = currentNode->pathDistance + 1;

                    priorityQueue.push(newNode);
                }
            }
        }

        lastDistanceToGoal = currentNode->distanceToGoal;
        endPoint = currentNode->point;
        if (currentNode->distanceToGoal > initialDistanceToGoal) {
            //No progress
            return;
        }

        int max = 20;
        movedCloserToGoal = true;
        while (currentNode->parent != nullptr && max > 0) {
            outDirections.push_back(currentNode->direction);
            currentNode = currentNode->parent;
            max--;
        }

        std::reverse(outDirections.begin(), outDirections.end());
    }
};
//======================================================================================================================

//======================================================================================================================
//  MOVE EVALUATION IMPLEMENTATION
//======================================================================================================================
class IMoveEvaluator {
public:
    virtual double Evaluate(const GameState &gameState, const string &move, int playerId) const = 0;
};

class RandomMoveEvaluator : public IMoveEvaluator {
public:
    double Evaluate(const GameState &gameState, const string &move, const int playerId) const override {
        return rand() % 1000;
    }
};

class LengthBasedMoveEvaluator : public IMoveEvaluator {
public:
    double Evaluate(const GameState &gameState, const string &move, const int playerId) const override {
        return move.length();
    }
};

class ExecutionMoveEvaluator : public IMoveEvaluator {
private:
    const double completeQuestCoef = 1;
    const double exitPointsCoef = 0.7;
    const double invalidMovementCoef = -1;
    const double proximityToQuestsCoef = 0.8;
    const double minimumMovementsCoef = 0.1;

public:
    double Evaluate(const GameState &gameState, const string &move, const int playerId) const override {
        string map = gameState.map;
        Point playerPosition = gameState.playerPositions.at(playerId);
        string playerTile = Map::GetTileAt(map, playerPosition);
        vector<quest> playerQuests = gameState.quests.at(playerId);

        double score = 0;
        for (int i = 0; i < move.size(); ++i) {
            Direction direction = Direction::UP;

            switch (move[i]) {
                case 'R' :
                    direction = Direction::RIGHT;
                    break;
                case 'D':
                    direction = Direction::DOWN;
                    break;
                case 'L' :
                    direction = Direction::LEFT;
                    break;
            }

            const Point newPosition = Point::GetPointInDirectionOf(playerPosition, direction);
            string newTile = Map::GetTileAt(map, newPosition);
            MapTile newMapTile(newTile, newPosition);

            Direction connectedDirection;
            if (!MapTile(playerTile, playerPosition).IsConnectedToTile(newMapTile, connectedDirection)) {
                //Invalid
                score += 1 * invalidMovementCoef;
            }

            for (quest aQuest: playerQuests) {
                Item item = gameState.items.at(aQuest);
                if (item.position == newPosition) {
                    score += 1 * completeQuestCoef;
                }
            }

            playerPosition = newPosition;
            playerTile = newTile;
        }

        score += std::count(playerTile.begin(), playerTile.end(), '1') * exitPointsCoef;
        score += 1 / move.size() * minimumMovementsCoef;

        for (quest aQuest: playerQuests) {
            Item item = gameState.items.at(aQuest);
            score += 1 / (Distance::GetDistance(item.position, playerPosition) + 1) * proximityToQuestsCoef;
        }

        return score;
    }
};

//======================================================================================================================
//======================================================================================================================
//  Genetic Evolution Implementation
//======================================================================================================================
class GASolver {
public:
    struct GAParams {
        int populationSize{15};
        int numberOfIterations{15};
        double mutationProbability{0.5};
    };

private:
    typedef string Gene;
    const GameState &gameState;
    IMoveEvaluator *moveEvaluator;
    const int playerId;
    const string directions = "URDL";
    GAParams params;

public:
    GASolver(IMoveEvaluator *moveEvaluator, const GameState &gameState, const int playerId) :
            moveEvaluator(moveEvaluator), gameState(gameState), playerId(playerId) {}

    Gene Solve() {
        cerr << "-----------------------------------------------------" << endl;
        cerr << "Solving with Genetic Evolution..." << endl;
        cerr << "Generating initial population..." << endl;

        vector<Gene> population = GenerateSample(params.populationSize);
        cerr << "Generated population. Size: " << population.size() << endl;
        for (int i = 0; i < params.numberOfIterations; ++i) {
            population = Evolve(population);
        }

        return population[0];
    }

private:
    vector<Gene> GenerateSample(int sampleSize) {
        vector<Gene> output;

        for (int i = 0; i < sampleSize; ++i) {
            output.push_back(GenerateGene());
        }

        return output;
    }

    Gene GenerateGene() {
        int length = (rand() % 20) + 1;
        Gene output;
        for (int i = 0; i < length; ++i) {
            output += directions[rand() % directions.size()];
        }

        return output;
    }

    vector<Gene> Evolve(vector<Gene> &currentGeneration) {
        //Get fittest individuals
        sort(currentGeneration.begin(), currentGeneration.end(), [this](Gene a, Gene b) {
            return GetFitness(a) > GetFitness(b);
        });

        //Ditch unfit individuals
        currentGeneration.resize(currentGeneration.size() / 2);

        //Get new generation
        vector<Gene> newGeneration;
        for (int i = 0; i < currentGeneration.size(); ++i) {
            Gene parent = currentGeneration[i];

            if (ShouldMutate()) {
                newGeneration.push_back(Mutate(parent));
            } else {
                int mateIndex = rand() % currentGeneration.size();
                // Avoid crossing same individuals
                while (mateIndex == i) {
                    mateIndex = rand() % currentGeneration.size();
                }
                Gene mate = currentGeneration[mateIndex];
                newGeneration.push_back(Crossover(parent, mate));
            }
        }

        //Add new generation to population
        currentGeneration.insert(std::end(currentGeneration), std::begin(newGeneration), std::end(newGeneration));
        return currentGeneration;
    }

    double GetFitness(Gene individual) {
        return moveEvaluator->Evaluate(gameState, individual, playerId);
    }

    bool ShouldMutate() {
        return ((double) rand() / (RAND_MAX)) < params.mutationProbability;
    }

    Gene Mutate(Gene individual) {
        int mutationIndex = rand() % individual.size();

        char replacementChar = directions[rand() % (directions.size() - 1)];

        //Ensure replacement is different from current
        while (individual[mutationIndex] == replacementChar) {
            replacementChar = directions[rand() % (directions.size() - 1)];
        }

        string replacement = string(1, replacementChar);

        Gene output = individual.replace(mutationIndex, 1, replacement);
        return output;

    }

    Gene Crossover(Gene parentA, Gene parentB) {
        const int &size = parentB.size() / 2;
        Gene bReplacement = parentB.substr(size, size + 1);

        parentA.resize(parentA.size() / 2);
        parentA += bReplacement;

        return parentA;
    }
};

//======================================================================================================================

//======================================================================================================================
//  MOVEMENT SOLUTION IMPLEMENTATION
//======================================================================================================================
class IMover {
public:
    virtual vector<Direction> GetDirections(const GameState &gameState, int playerId) const = 0;
};

class AStarMover : public IMover {
public:
    vector<Direction> GetDirections(const GameState &gameState, const int playerId) const override {
        vector<quest> quests = gameState.quests.at(playerId);
        Point playerPosition = gameState.playerPositions.at(playerId);
        AStar aStar;

        vector<Direction> bestDirection;
        int closestDistanceToGoal = 49;

        AStarResult output;
        //TODO: Try passing in a list of bonus goals to astar and see if resulting path crossed any of those.
        for (const quest &theQuest: quests) {
            output = aStar.DoSearch(gameState.map, playerPosition, gameState.items.at(theQuest).position);

            if (output.reachedGoal) {
                return output.directions;
            }

            if (output.finalDistanceToGoal < closestDistanceToGoal) {
                closestDistanceToGoal = output.finalDistanceToGoal;
                bestDirection = output.directions;
            }
        }

        return bestDirection;
    }
};

class QuestBasedMover : public IMover {
    vector<Direction> GetDirections(const GameState &gameState, int playerId) const override {
        map<quest, bool> explored;
        vector<quest> quests = gameState.quests.at(playerId);

        vector<Direction> currentPath;
        AStar aStar;
        AStarResult aStarResult;
        Point lastPoint = gameState.playerPositions.at(playerId);
        int i = 0;
        vector<Direction> bestEndPath;
        int closestDistanceToGoal = 49;

        bool isRerunning = false;
        while (currentPath.size() < 20 && explored.size() != quests.size()) {
            if (i >= quests.size() && !isRerunning) {
                i = 0;
                isRerunning = true;
            } else if (i >= quests.size()) {
                break;
            }

            quest targetQuest = quests[i];
            if (explored.find(targetQuest) != explored.end()) {
                i++;
                continue;
            }

            aStarResult = aStar.DoSearch(gameState.map, lastPoint, gameState.items.at(targetQuest).position);
            if (aStarResult.reachedGoal && (currentPath.size() + aStarResult.directions.size()) <= 20) {
                currentPath.insert(std::end(currentPath), std::begin(aStarResult.directions),
                                   std::end(aStarResult.directions));
                explored[targetQuest] = true;
                lastPoint = aStarResult.endPoint;
            } else if (isRerunning && aStarResult.finalDistanceToGoal < closestDistanceToGoal) {
                closestDistanceToGoal = aStarResult.finalDistanceToGoal;
                bestEndPath = aStarResult.directions;
                explored[targetQuest] = true;
            }
            i++;
        }

        int spaceLeft = 20 - currentPath.size();
        if (spaceLeft < bestEndPath.size()) {
            bestEndPath.resize(spaceLeft);
        }

        currentPath.insert(std::end(currentPath), std::begin(bestEndPath), std::end(bestEndPath));
        return currentPath;
    }
};

class GAMover : public IMover {
    vector<Direction> GetDirections(const GameState &gameState, const int playerId) const override {
        GASolver solver(new ExecutionMoveEvaluator(), gameState, playerId);
        string directions = solver.Solve();
        vector<Direction> output;
        for (int i = directions.size() - 1; i >= 0; --i) {
            switch (directions[i]) {
                case 'U':
                    output.push_back(Direction::UP);
                    break;
                case 'R':
                    output.push_back(Direction::RIGHT);
                    break;
                case 'D':
                    output.push_back(Direction::DOWN);
                    break;
                case 'L':
                    output.push_back(Direction::LEFT);
                    break;
            }
        }

        return output;
    }
};

//======================================================================================================================
//  GAME EVALUATION IMPLEMENTATION
//======================================================================================================================
class IEvaluator {
public:
    virtual double Evaluate(const GameState &gameState, int playerId) const = 0;
};

class GreedyEvaluator : public IEvaluator {
private:
    const double distanceToGoalCoef = 0.35;
    const double existingPathToGoalCoef = 1;
    const double existingPathCloseToGoalCoef = 0.8;
    const double oppDistanceToGoalCoef = -0.3;
    const double oppExistingPathToGoalCoef = -0.9;
    const double oppExistingPathCloseToGoalCoef = -0.85;

public:
    double Evaluate(const GameState &gameState, const int playerId) const override {
        int oppId = playerId == 1 ? 0 : 1;

        Point playerPosition = gameState.playerPositions.at(playerId);
        Point oppPosition = gameState.playerPositions.at(oppId);

        AStar aStar;
        AStarResult aStarResult;
        double score = 0;
        double distanceScore;
        double pathToGoalScore;
        double closerToGoalScore;

        for (const quest &playerQuest: gameState.quests.at(playerId)) {
            Point playerQuestPosition = gameState.items.at(playerQuest).position;
            aStarResult = aStar.DoSearch(gameState.map, playerPosition, playerQuestPosition);
            distanceScore =
                    1.0 / (Distance::GetDistance(playerPosition, playerQuestPosition) + 0.000000001) *
                    distanceToGoalCoef;
            pathToGoalScore = (double) aStarResult.reachedGoal * existingPathToGoalCoef;
            closerToGoalScore = 1.0 / (aStarResult.finalDistanceToGoal + 0.000000001) * existingPathCloseToGoalCoef;
            score += distanceScore;
            score += pathToGoalScore;
            score += closerToGoalScore;
        }

        for (const quest &oppQuest: gameState.quests.at(playerId)) {
            Point oppQuestPosition = gameState.items.at(oppQuest).position;
            aStarResult = aStar.DoSearch(gameState.map, oppPosition, oppQuestPosition);

            distanceScore =
                    1.0 / (Distance::GetDistance(oppPosition, oppQuestPosition) + 0.000000001) * oppDistanceToGoalCoef;
            pathToGoalScore = (double) aStarResult.reachedGoal * oppExistingPathToGoalCoef;
            closerToGoalScore = 1.0 / (aStarResult.finalDistanceToGoal + 0.000000001) * oppExistingPathCloseToGoalCoef;
            score += distanceScore;
            score += pathToGoalScore;
            score += closerToGoalScore;
        }

        return score;
    }
};

//======================================================================================================================
//  PUSH SOLUTION IMPLEMENTATION
//======================================================================================================================
class IPusher {
public:
    virtual pair<int, Direction> Push(const GameState &gameState, int playerId) const = 0;
};

class BruteForcePusher : public IPusher {
public:
    pair<int, Direction> Push(const GameState &gameState, const int playerId) const override {
        IEvaluator *evaluator = new GreedyEvaluator();
        pair<int, Direction> bestSolution;
        double bestScore = 0;

        for (int rowColId = 0; rowColId < 7; ++rowColId) {
            for (auto direction: Dir::Directions) {
                GameState pushedState = gameState.Push(rowColId, direction, playerId);

                double score = evaluator->Evaluate(pushedState, playerId);

                if (score > bestScore) {
                    bestScore = score;
                    bestSolution = {rowColId, direction};
                }
            }
        }

        return bestSolution;
    }
};
//======================================================================================================================

class Solver {
public:
    static void DoPush(const GameState &gameState, const IPusher *pusher) {
        pair<int, Direction> result = pusher->Push(gameState, ACTOR_ID);
        cout << "PUSH " << result.first << " " << Dir::GetAsString(result.second) << endl;
    }

    static void DoMove(const GameState &gameState, const IMover *mover) {
        vector<Direction> directions = mover->GetDirections(gameState, ACTOR_ID);

        if (directions.empty()) {
            cout << "PASS" << endl;
            return;
        }

        cout << "MOVE";
        for (int i = 0; i < directions.size(); ++i) {
            Direction direction = directions[i];
            cout << " " << Dir::GetAsString(direction);
        }

        cout << endl;
    }
};

/**
 * Help the Christmas elves fetch presents in a magical labyrinth!
 **/
int main() {

//    IMover *mover = new AStarMover();
//    IMover *mover = new GAMover();
    IMover *mover = new QuestBasedMover();
    IPusher *pusher = new BruteForcePusher();

    // game loop
    while (true) {
        int turnType;
        cin >> turnType;
        cin.ignore();

        GameState gameState;
        //------------------- READ MAP ------------------------------------------------------------
        string mapStr;
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 7; j++) {
                string tile;
                cin >> tile;
                cin.ignore();
                mapStr += tile;
            }
        }

        gameState.map = mapStr;

        //----------------- READ PLAYER INFO -----------------------------------------------------
        Point actorPosition;
        Point villainPosition;
        string actorTile;
        string villainTile;
        for (int playerId = 0; playerId < 2; playerId++) {
            int numPlayerCards; // the total number of quests for a player (hidden and revealed)
            int playerX;
            int playerY;
            string playerTile;
            cin >> numPlayerCards >> playerX >> playerY >> playerTile;
            cin.ignore();

            if (playerId == ACTOR_ID) {
                actorPosition = Point(playerX, playerY);
                actorTile = playerTile;
                cerr << "Actor is at: " << actorPosition.toString() << "(" << actorTile << ")" << endl;
            } else {
                villainPosition = Point(playerX, playerY);
                villainTile = playerTile;
                cerr << "Villain is at: " << villainPosition.toString() << "(" << villainTile << ")" << endl;
            }
        }

        gameState.playerPositions[ACTOR_ID] = actorPosition;
        gameState.playerPositions[VILLAIN_ID] = villainPosition;
        gameState.playerTiles[ACTOR_ID] = actorTile;
        gameState.playerTiles[VILLAIN_ID] = villainTile;

        //------------------ READ ITEMS -----------------------------------------------------
        int numItems; // the total number of items available on board and on player tiles
        cin >> numItems;
        cin.ignore();
        std::map<quest, Item> items;
        for (int i = 0; i < numItems; i++) {
            string itemName;
            int itemX;
            int itemY;
            int itemPlayerId;
            cin >> itemName >> itemX >> itemY >> itemPlayerId;
            const Item &item = Item(itemName, itemX, itemY, itemPlayerId);
            items[quest(itemName, itemPlayerId)] = item;
            cin.ignore();
        }

        gameState.items = items;

        //------------------- READ QUESTS ----------------------------------------------------
        int numQuests; // the total number of revealed quests for both players
        cin >> numQuests;
        cin.ignore();

        map<int, vector<quest>> quests;

        for (int i = 0; i < numQuests; i++) {
            string questItemName;
            int questPlayerId;
            cin >> questItemName >> questPlayerId;
            cin.ignore();

            quest theQuest = quest(questItemName, questPlayerId);
            if (quests.find(questPlayerId) == quests.end()) {
                quests[questPlayerId] = vector<quest>();
            }

            quests[questPlayerId].push_back(theQuest);
        }

        gameState.quests = quests;

        // ------------------ GET SOLUTION ----------------------------------------------------

        if (turnType == TURN_TYPE_PUSH) {
            Solver::DoPush(gameState, pusher);
        } else {
            Solver::DoMove(gameState, mover);
        }
    }
}
