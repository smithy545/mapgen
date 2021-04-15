//
// Created by Philip Smith on 9/16/2020.
//

#ifndef MAPGEN_FORTUNEALGORITHM_H
#define MAPGEN_FORTUNEALGORITHM_H

#include <glm/glm.hpp>
#include <unordered_map>
#include <math_util.h>
#include <string>
#include <vector>

#include "Diagram.h"


class FortuneAlgorithm {
public:
    struct Event;

    struct HalfEdge {
        glm::vec2 origin;
        glm::vec2 destination;
        HalfEdge *prev, *next;
    };

    struct Arc {
        glm::vec2 focus;
        HalfEdge *prev_edge, *next_edge;
        Arc *parent;
        Arc *left, *right;
        Arc *prev, *next;
        enum {
            BLACK,
            RED
        } color; // for red/black tree balancing
        Event *event{nullptr}; // for faster removal
    };

    struct Cell {
        HalfEdge *handle{nullptr};
    };

    struct Event {
        enum {
            SITE,
            CIRCLE
        } type;
        double y;
        glm::vec2 position;
        Arc *arc;
        Event *prev{nullptr}, *next{nullptr};
    };

    class EventQueue {
    public:
        void insert(Event *event);

        void remove(Arc *arc);

        bool empty();

        Event *pop();

        ~EventQueue();
    private:
        Event *head{nullptr}, *tail{nullptr};
    };

    class Beachline {
    public:
        Arc* create_arc(glm::vec2 focus);

        void delete_arc_tree(Arc* tree_root);

        void set_root(glm::vec2 site);

        Arc get_seam();

        void insert_before(Arc *parent, Arc *child);

        void insert_after(Arc *parent, Arc *child);

        void insert_fixup(Arc *x);

        void delete_from(Arc *x);

        void delete_fixup(Arc *x);

        void replace(Arc *old, Arc *replacement);

        Arc *minimum(Arc *x);

        Arc *search(double x, double directrix_y);

        void left_rotate(Arc *x);

        void right_rotate(Arc *x);

        void transplant(Arc *old, Arc *replacement);

        bool is_nil(Arc *x);

        void print();

        ~Beachline();
    private:
        Arc *root{nullptr};
        Arc sentinel{
                glm::vec2(),
                nullptr,
                nullptr,
                &sentinel,
                &sentinel,
                &sentinel,
                &sentinel,
                &sentinel,
                Arc::BLACK
        };
        Arc seam{
                glm::vec2(),
                nullptr,
                nullptr,
                &sentinel,
                &sentinel,
                &sentinel,
                &sentinel,
                &sentinel,
                Arc::BLACK
        };
        double epsilon{1.0}; // ~zero for use in double comparisons
    };

    explicit FortuneAlgorithm(const std::vector<glm::vec2>& sites);

    ~FortuneAlgorithm();

    Diagram construct(bool validate_diagram = true);

    Diagram construct(double width, double height, bool validate_diagram = true, double x_offset = 0.0, double y_offset = 0.0);

    static bool voroni_sort(glm::vec2 v1, glm::vec2 v2);

    static bool check_for_edge_intersections(std::vector<Diagram::Edge> edges);
private:
    std::unordered_map<std::string, Cell> cells;
    Beachline beachline;
    EventQueue events;
    double sweep_y;

    void generate_cell(glm::vec2 site, HalfEdge* edge_handle);

    void split_arc(Arc *arc, glm::vec2 site);

    void handle_site(Event *event);

    void handle_circle(Event *event);

    bool add_circle_event(Arc *arc);

    static HalfEdge* create_half_edge(glm::vec2 origin);

    // math
    glm::vec2 compute_edge_origin(glm::vec2 parent, glm::vec2 child);

    static glm::vec2 compute_center(glm::vec2 left, glm::vec2 middle, glm::vec2 right);
};

#endif //MAPGEN_FORTUNEALGORITHM_H
