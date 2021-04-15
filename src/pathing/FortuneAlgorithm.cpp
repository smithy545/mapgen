//
// Created by Philip Smith on 9/16/2020.
//

#include "FortuneAlgorithm.h"

#include <cfloat>
#include <fmt/format.h>
#include <iostream>
#include <limits>
#include <unordered_set>


FortuneAlgorithm::FortuneAlgorithm(const std::vector<glm::vec2>& sites) {
    assert(!sites.empty());
    for (auto site: sites)
        events.insert(new Event{Event::SITE, site.y, site});
}

FortuneAlgorithm::~FortuneAlgorithm() {
    std::unordered_set<HalfEdge*> deleted_edges;
    for(auto [site, cell]: cells) {
        auto original_handle = cell.handle;
        auto handle = original_handle;
        do {
            auto prev = handle;
            handle = handle->next;
            prev->next = nullptr;
            deleted_edges.insert(prev);
            if(!deleted_edges.contains(prev))
                delete prev;
        } while(handle != nullptr && handle != original_handle);
    }
}

Diagram FortuneAlgorithm::construct(bool validate_diagram) {
    auto initial_event = events.pop();
    beachline.set_root(initial_event->position);
    while (!events.empty()) {
        auto event = events.pop();
        sweep_y = event->y;
        switch (event->type) {
            case Event::SITE:
                handle_site(event);
                break;
            case Event::CIRCLE:
                handle_circle(event);
        }
        delete event;
    }

    // ignore unfinished half edges
    Diagram voroni;
    std::unordered_set<std::string> edge_set;
    for(auto [site, cell]: cells) {
        auto original_handle = cell.handle;
        auto handle = original_handle;
        do {
            if(!edge_set.contains(edge_key(handle->destination, handle->origin))) {
                edge_set.insert(edge_key(handle->origin, handle->destination));
                voroni.edges.emplace_back(handle->origin, handle->destination);
            }
            handle = handle->next;
        } while(handle != nullptr && handle != original_handle);
    }

    if(validate_diagram && check_for_edge_intersections(voroni.edges))
        throw "Edge intersection found. Voroni diagram invalid.";

    return voroni;
}

Diagram FortuneAlgorithm::construct(double width, double height, bool validate_diagram, double x_offset, double y_offset) {
    auto initial_event = events.pop();
    beachline.set_root(initial_event->position);
    while (!events.empty()) {
        auto event = events.pop();
        sweep_y = event->y;
        switch (event->type) {
            case Event::SITE:
                handle_site(event);
                break;
            case Event::CIRCLE:
                handle_circle(event);
        }
        delete event;
    }

    // cut off half edges outside boundary and extend edges within boundary to edge if necessary
    Diagram voroni;
    std::unordered_set<std::string> edge_set;
    utils::math::rect bounds{x_offset, y_offset, width, height};
    for(auto [site, cell]: cells) {
        auto original_handle = cell.handle;
        auto handle = original_handle;
        do {
            if(!edge_set.contains(edge_key(handle->destination, handle->origin)) &&
            (utils::math::in_rect(handle->origin, bounds)
            || utils::math::in_rect(handle->destination, bounds)
            || utils::math::do_intersect(handle->origin, handle->destination, {x_offset, y_offset}, {x_offset+width, y_offset})
            || utils::math::do_intersect(handle->origin, handle->destination, {x_offset+width, y_offset}, {x_offset+width, y_offset+height})
            || utils::math::do_intersect(handle->origin, handle->destination, {x_offset+width, y_offset+height}, {x_offset, y_offset+height})
            || utils::math::do_intersect(handle->origin, handle->destination, {x_offset, y_offset+height}, {x_offset, y_offset}))) {
                edge_set.insert(edge_key(handle->origin, handle->destination));
                voroni.edges.emplace_back(handle->origin, handle->destination);
            }
            handle = handle->next;
        } while(handle != nullptr && handle != original_handle);
    }
    voroni.edges.emplace_back(glm::vec2{x_offset+1, y_offset+1}, glm::vec2{x_offset+width, y_offset+1});
    voroni.edges.emplace_back(glm::vec2{x_offset+width, y_offset+1}, glm::vec2{x_offset+width, y_offset+height-1});
    voroni.edges.emplace_back(glm::vec2{x_offset+width, y_offset+height-1}, glm::vec2{x_offset+1, y_offset+height-1});
    voroni.edges.emplace_back(glm::vec2{x_offset+1, y_offset+height-1}, glm::vec2{x_offset+1, y_offset+1});

    if(validate_diagram && check_for_edge_intersections(voroni.edges))
        throw "Edge intersection found. Voroni diagram invalid.";

    return voroni;
}

// event handlers
void FortuneAlgorithm::handle_site(Event *event) {
    auto site = event->position;
    auto parent = beachline.search(site.x, sweep_y);
    split_arc(parent, site);
}

void FortuneAlgorithm::handle_circle(Event *event) {
    auto position = event->position;
    auto arc = event->arc;
    auto left = arc->prev;
    auto right = arc->next;
    auto left_edge = arc->prev_edge;
    auto right_edge = arc->next_edge;

    // delete squeezed arc and remove events
    events.remove(left);
    events.remove(right);
    beachline.delete_from(arc);
    delete arc;

    // update edge positions
    left_edge->destination = position;
    right_edge->origin = position;
    left->next_edge->origin = position;
    right->prev_edge->destination = position;

    // update edge relations
    left_edge->next = right_edge;
    right_edge->prev = left_edge;
    // left
    left->next_edge->prev = create_half_edge(position);
    left->next_edge->prev->next = left->next_edge;
    left->next_edge = left->next_edge->prev;
    // right
    right->prev_edge->next = create_half_edge(position);
    right->prev_edge->next->prev = right->prev_edge;
    right->prev_edge = right->prev_edge->next;

    // new circle events
    if (!beachline.is_nil(left->prev))
        add_circle_event(left);
    if (!beachline.is_nil(right->next))
        add_circle_event(right);
}

void FortuneAlgorithm::split_arc(Arc *arc, glm::vec2 site) {
    auto child = beachline.create_arc(site);
    auto seam = beachline.get_seam();
    Arc *left, *right;
    if(!beachline.is_nil(seam.next) && !beachline.is_nil(seam.prev)) {
        // if new arc lands on breakpoint between existing arcs treat as a circle event then create edge for each parent
        left = seam.prev;
        right = seam.next;
        auto edge_origin = compute_center(left->focus, site, right->focus);
        left->next_edge->origin = edge_origin;
        right->prev_edge->destination = edge_origin;

        events.remove(left);
        events.remove(right);
        beachline.insert_before(right, child);

        // create new pairs
        left->next_edge->prev = create_half_edge(edge_origin);
        left->next_edge->prev->next = left->next_edge;
        left->next_edge = left->next_edge->prev;
        child->prev_edge = create_half_edge(edge_origin);

        right->prev_edge->next = create_half_edge(edge_origin);
        right->prev_edge->next->prev = right->prev_edge;
        right->prev_edge = right->prev_edge->next;
        child->next_edge = create_half_edge(edge_origin);

        child->prev_edge->prev = child->next_edge;
        child->next_edge->next = child->prev_edge;
        generate_cell(site, child->prev_edge);
    } else if(arc->focus.y == site.y) { // special case handling for better precision
        // always on right side due to site sorting
        left = arc;
        right = arc->next;
        child->prev = arc;
        child->next = arc->next;
        auto x = 0.5f * (site.x + arc->focus.x);
        auto a = (x - arc->focus.x);
        auto edge_origin = glm::vec2(x, sweep_y - std::sqrt(2.f) * a);
        child->prev_edge = create_half_edge(edge_origin);
        if(!beachline.is_nil(child->next)) {
            child->next->prev = child;
            child->next->prev_edge = create_half_edge(compute_edge_origin(child->next->focus, child->focus));
            child->next_edge = child->prev_edge;
        }
        arc->next_edge = create_half_edge(edge_origin);
        events.remove(arc);
        beachline.insert_after(arc, child);
    } else { // split parent into two neighboring arcs, insert child between, create new edge shared by parent and child
        left = beachline.create_arc(arc->focus);
        right = beachline.create_arc(arc->focus);
        right->next_edge = arc->next_edge;
        left->prev_edge = arc->prev_edge;

        events.remove(arc);
        beachline.replace(arc, child);
        beachline.insert_before(child, left);
        beachline.insert_after(child, right);

        auto edge_origin = compute_edge_origin(arc->focus, site);
        delete arc;

        // create child-side edge
        child->prev_edge = create_half_edge(edge_origin);
        child->next_edge = child->prev_edge;
        // create parent-side edge
        right->prev_edge = create_half_edge(edge_origin);
        left->next_edge = right->prev_edge;
        // add cell if it doesn't exist
        generate_cell(site, child->prev_edge);
    }
    // add circle events
    if (!beachline.is_nil(left->prev))
        add_circle_event(left);
    if (!beachline.is_nil(right->next))
        add_circle_event(right);
}

bool FortuneAlgorithm::add_circle_event(Arc *arc) {
    auto left = arc->prev->focus;
    auto middle = arc->focus;
    auto right = arc->next->focus;
    auto collision = compute_center(left, middle, right);
    auto event_y = collision.y + glm::distance(left, collision);
    auto left_moving_right = left.y > middle.y;
    auto right_moving_left = middle.y > right.y;
    int cx = collision.x;
    int lx = left.x;
    int mx = middle.x;
    int rx = right.x;
    auto valid = ((left_moving_right && lx <= cx) || (!left_moving_right && mx > cx))
            && ((right_moving_left && mx < cx) || (!right_moving_left && rx >= cx));
    if (event_y >= sweep_y && valid) {
        arc->event = new Event{Event::CIRCLE, event_y, collision, arc};
        events.insert(arc->event);
        return true;
    }
    return false;
}

// helper methods
FortuneAlgorithm::HalfEdge *FortuneAlgorithm::create_half_edge(glm::vec2 origin) {
    return new HalfEdge{
            origin,
            origin,
            nullptr,
            nullptr
    };
}

glm::vec2 FortuneAlgorithm::compute_edge_origin(glm::vec2 parent, glm::vec2 child) {
    auto origin = glm::vec2{child.x, utils::math::compute_parabola_y(parent, sweep_y, child.x)};
    return origin;
}

glm::vec2 FortuneAlgorithm::compute_center(glm::vec2 left, glm::vec2 middle, glm::vec2 right) {
    // Reference: https://github.com/pvigier/FortuneAlgorithm/blob/master/src/FortuneAlgorithm.cpp: line 212
    auto v1 = left - middle;
    auto v2 = middle - right;
    // orthogonal to edges
    v1 = glm::vec2(-v1.y, v1.x);
    v2 = glm::vec2(-v2.y, v2.x);
    // half edge
    auto delta = 0.5f * (right - left);
    float t = (delta.x * v2.y - delta.y * v2.x) / (v1.x * v2.y - v1.y * v2.x);
    return 0.5f * (left + middle) + t * v1;
}

std::string FortuneAlgorithm::site_key(glm::vec2 site, bool plain) {
    if(plain)
        return fmt::format("{},{}", (int)site.x, (int)site.y);
    return fmt::format("glm::vec2({},{})", (int)site.x, (int)site.y);
}

std::string FortuneAlgorithm::edge_key(glm::vec2 origin, glm::vec2 destination) {
    return fmt::format("{}, {} | {}, {}", origin.x, origin.y, destination.x, destination.y);
}

bool FortuneAlgorithm::check_for_edge_intersections(std::vector<FortuneAlgorithm::Edge> edges) {
    for(int i = 0; i < edges.size() - 1; i++) {
        for(int j = i + 1; j < edges.size(); j++) {
            if(edges[i].first.y >= 0 && edges[i].second.y >= 0 && edges[j].first.y >= 0 && edges[j].second.y >= 0
               && edges[i].first.x >= 0 && edges[i].second.x >= 0 && edges[j].first.x >= 0 && edges[j].second.x >= 0
               && glm::length(edges[i].first - edges[j].first) > 2.0
               && glm::length(edges[i].first - edges[j].second) > 2.0
               && glm::length(edges[i].second - edges[j].first) > 2.0
               && glm::length(edges[i].second - edges[j].second) > 2.0
               && utils::math::do_intersect(edges[i].first, edges[i].second, edges[j].first, edges[j].second)) {
                std::cout << "Error detected between the following edges" << std::endl;
                std::cout << edge_key(edges[i].first, edges[i].second) << std::endl;
                std::cout << edge_key(edges[j].first, edges[j].second) << std::endl;
                return true;
            }
        }
    }
    return false;
}

bool FortuneAlgorithm::voroni_sort(glm::vec2 v1, glm::vec2 v2) {
    return v1.y < v2.y || (v1.y == v2.y && v1.x < v2.x);
}

void FortuneAlgorithm::generate_cell(glm::vec2 site, FortuneAlgorithm::HalfEdge *edge_handle) {
    auto sk = site_key(site);
    if(!cells.contains(sk))
        cells[sk] = Cell{edge_handle};
}

// Beachline methods
FortuneAlgorithm::Arc *FortuneAlgorithm::Beachline::search(double x, double directrix_y) {
    auto itr = root;
    while (!is_nil(itr)) {
        double left_collision = -std::numeric_limits<double>::infinity();
        double right_collision = std::numeric_limits<double>::infinity();
        if (!is_nil(itr->prev)) {
            if(itr->prev->focus.y == directrix_y && itr->focus.y == directrix_y)
                left_collision = 0.5*(itr->focus.x + itr->prev->focus.x);
            else
                left_collision = utils::math::compute_parabolic_collision_x(itr->prev->focus, itr->focus, directrix_y);
        }
        if (!is_nil(itr->next)) {
            if(itr->next->focus.y == directrix_y && itr->focus.y == directrix_y)
                right_collision = 0.5*(itr->focus.x + itr->next->focus.x);
            else
                right_collision = utils::math::compute_parabolic_collision_x(itr->focus, itr->next->focus, directrix_y);
        }
        if (x < left_collision)
            itr = itr->left;
        else if (x > right_collision)
            itr = itr->right;
        else {
            // check if landing on the "seam" between two arcs
            if(std::abs(left_collision - x) < epsilon) {
                seam.prev = itr->prev;
                seam.next = itr;
            } else if(std::abs(right_collision - x) < epsilon) {
                seam.prev = itr;
                seam.next = itr->next;
            } else {
                seam.prev = &sentinel;
                seam.next = &sentinel;
            }
            return itr;
        }
    }
    return itr;
}

void FortuneAlgorithm::Beachline::insert_before(Arc *parent, Arc *child) {
    // Reference: https://github.com/pvigier/FortuneAlgorithm/blob/master/src/Beachline.cpp: line 87
    if (is_nil(parent->left)) {
        parent->left = child;
        child->parent = parent;
    } else {
        parent->prev->right = child;
        child->parent = parent->prev;
    }
    child->prev = parent->prev;
    if (!is_nil(child->prev))
        child->prev->next = child;
    child->next = parent;
    parent->prev = child;
    insert_fixup(child);
}

void FortuneAlgorithm::Beachline::insert_after(Arc *parent, Arc *child) {
    // Reference: https://github.com/pvigier/FortuneAlgorithm/blob/master/src/Beachline.cpp: line 110
    if (is_nil(parent->right)) {
        parent->right = child;
        child->parent = parent;
    } else {
        parent->next->left = child;
        child->parent = parent->next;
    }
    child->next = parent->next;
    if (!is_nil(child->next))
        child->next->prev = child;
    child->prev = parent;
    parent->next = child;
    insert_fixup(child);
}

void FortuneAlgorithm::Beachline::insert_fixup(Arc *x) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 316 RB-INSERT-FIXUP
    while (x->parent->color == Arc::RED) {
        if (x->parent == x->parent->parent->left) {
            auto y = x->parent->parent->right;
            if (y->color == Arc::RED) {
                x->parent->color = Arc::BLACK;
                y->color = Arc::BLACK;
                x->parent->parent->color = Arc::RED;
                x = x->parent->parent;
            } else {
                if (x == x->parent->right) {
                    x = x->parent;
                    left_rotate(x);
                }
                x->parent->color = Arc::BLACK;
                x->parent->parent->color = Arc::RED;
                right_rotate(x->parent->parent);
            }
        } else {
            auto y = x->parent->parent->left;
            if (y->color == Arc::RED) {
                x->parent->color = Arc::BLACK;
                y->color = Arc::BLACK;
                x->parent->parent->color = Arc::RED;
                x = x->parent->parent;
            } else {
                if (x == x->parent->left) {
                    x = x->parent;
                    right_rotate(x);
                }
                x->parent->color = Arc::BLACK;
                x->parent->parent->color = Arc::RED;
                left_rotate(x->parent->parent);
            }
        }
    }
    root->color = Arc::BLACK;
}

void FortuneAlgorithm::Beachline::delete_from(Arc *z) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 324 RB-DELETE
    auto y = z;
    auto y_original_color = y->color;
    Arc *x;
    if (is_nil(z->left)) {
        x = z->right;
        transplant(z, z->right);
    } else if (is_nil(z->right)) {
        x = z->left;
        transplant(z, z->left);
    } else {
        y = minimum(z->right);
        y_original_color = y->color;
        x = y->right;
        if (y->parent == z)
            x->parent = y;
        else {
            transplant(y, y->right);
            y->right = z->right;
            y->right->parent = y;
        }
        transplant(z, y);
        y->left = z->left;
        y->left->parent = y;
        y->color = z->color;
    }
    if (y_original_color == Arc::BLACK)
        delete_fixup(x);
    if (!is_nil(z->prev))
        z->prev->next = z->next;
    if (!is_nil(z->next))
        z->next->prev = z->prev;
}

void FortuneAlgorithm::Beachline::delete_fixup(Arc *x) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 326 RB-DELETE-FIXUP
    while (x != root && x->color == Arc::BLACK) {
        if (x == x->parent->left) {
            auto w = x->parent->right;
            if (w->color == Arc::RED) {
                w->color = Arc::BLACK;
                x->parent->color = Arc::RED;
                left_rotate(x->parent);
                w = x->parent->right;
            }
            if (w->left->color == Arc::BLACK && w->right->color == Arc::BLACK) {
                w->color = Arc::RED;
                x = x->parent;
            } else {
                if (w->right->color == Arc::BLACK) {
                    w->left->color = Arc::BLACK;
                    w->color = Arc::RED;
                    right_rotate(w);
                    w = x->parent->right;
                }
                w->color = x->parent->color;
                x->parent->color = Arc::BLACK;
                w->right->color = Arc::BLACK;
                left_rotate(x->parent);
                x = root;
            }
        } else {
            auto w = x->parent->left;
            if (w->color == Arc::RED) {
                w->color = Arc::BLACK;
                x->parent->color = Arc::RED;
                right_rotate(x->parent);
                w = x->parent->left;
            }
            if (w->left->color == Arc::BLACK && w->right->color == Arc::BLACK) {
                w->color = Arc::RED;
                x = x->parent;
            } else {
                if (w->left->color == Arc::BLACK) {
                    w->right->color = Arc::BLACK;
                    w->color = Arc::RED;
                    left_rotate(w);
                    w = x->parent->left;
                }
                w->color = x->parent->color;
                x->parent->color = Arc::BLACK;
                w->left->color = Arc::BLACK;
                right_rotate(x->parent);
                x = root;
            }
        }
    }
    x->color = Arc::BLACK;
}

void FortuneAlgorithm::Beachline::replace(FortuneAlgorithm::Arc *old, FortuneAlgorithm::Arc *rep) {
    // Reference: https://github.com/pvigier/FortuneAlgorithm/blob/master/src/Beachline.cpp: line 133
    transplant(old, rep);
    rep->left = old->left;
    rep->right = old->right;
    if (!is_nil(rep->left))
        rep->left->parent = rep;
    if (!is_nil(rep->right))
        rep->right->parent = rep;
    rep->prev = old->prev;
    rep->next = old->next;
    if (!is_nil(rep->prev))
        rep->prev->next = rep;
    if (!is_nil(rep->next))
        rep->next->prev = rep;
    rep->color = old->color;
}

void FortuneAlgorithm::Beachline::left_rotate(Arc *x) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 312 LEFT-ROTATE
    auto y = x->right;
    x->right = y->left;
    if (!is_nil(y->left))
        y->left->parent = x;
    y->parent = x->parent;
    if (is_nil(x->parent))
        root = y;
    else if (x == x->parent->left)
        x->parent->left = y;
    else
        x->parent->right = y;
    y->left = x;
    x->parent = y;
}

void FortuneAlgorithm::Beachline::right_rotate(Arc *x) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 312 LEFT-ROTATE
    auto y = x->left;
    x->left = y->right;
    if (!is_nil(y->right))
        y->right->parent = x;
    y->parent = x->parent;
    if (is_nil(x->parent))
        root = y;
    else if (x == x->parent->right)
        x->parent->right = y;
    else
        x->parent->left = y;
    y->right = x;
    x->parent = y;
}

void FortuneAlgorithm::Beachline::transplant(Arc *old, Arc *rep) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 296 TRANSPLANT
    if (is_nil(old->parent))
        root = rep;
    else if (old == old->parent->left)
        old->parent->left = rep;
    else
        old->parent->right = rep;
    rep->parent = old->parent;
}

FortuneAlgorithm::Arc *FortuneAlgorithm::Beachline::minimum(Arc *x) {
    // Reference: Introduction to Algorithms by Cormen, Leiserson, Rivest, Stein: page 291 TREE-MINIMUM
    while (!is_nil(x->left))
        x = x->left;
    return x;
}

bool FortuneAlgorithm::Beachline::is_nil(FortuneAlgorithm::Arc *x) {
    return x == &sentinel;
}

FortuneAlgorithm::Arc *FortuneAlgorithm::Beachline::create_arc(glm::vec2 focus) {
    return new Arc{
            focus,
            nullptr,
            nullptr,
            &sentinel,
            &sentinel,
            &sentinel,
            &sentinel,
            &sentinel,
            Arc::RED
    };
}

void FortuneAlgorithm::Beachline::delete_arc_tree(FortuneAlgorithm::Arc *tree_root) {
    if(is_nil(tree_root))
        return;
    delete tree_root->prev_edge;
    delete tree_root->next_edge;
    delete_arc_tree(tree_root->left);
    delete_arc_tree(tree_root->right);
    delete tree_root;
}

void FortuneAlgorithm::Beachline::set_root(glm::vec2 site) {
    root = create_arc(site);
    root->color = Arc::BLACK;
}

FortuneAlgorithm::Arc FortuneAlgorithm::Beachline::get_seam() {
    return seam;
}

FortuneAlgorithm::Beachline::~Beachline() {
    delete_arc_tree(root);
}

void FortuneAlgorithm::Beachline::print() {
    auto itr = root;
    while(!is_nil(itr->prev))
        itr = itr->prev;
    while(!is_nil(itr->next)) {
        std::cout << site_key(itr->focus, true) << " - ";
        itr = itr->next;
    }
    std::cout << site_key(itr->focus, true) << std::endl;
}

// EventQueue methods
void FortuneAlgorithm::EventQueue::insert(FortuneAlgorithm::Event *event) {
    if(head == nullptr)
        head = event;
    else {
        auto itr = head;
        while (true) {
            if (itr->y > event->y ||
                (itr->y == event->y && (itr->type == Event::CIRCLE || event->position.x < itr->position.x))) {
                event->next = itr;
                event->prev = itr->prev;
                if (event->prev != nullptr)
                    event->prev->next = event;
                else
                    head = event;
                itr->prev = event;
                return;
            } else if (itr->next == nullptr) {
                itr->next = event;
                event->prev = itr;
                return;
            }
            itr = itr->next;
        }
    }
}

void FortuneAlgorithm::EventQueue::remove(FortuneAlgorithm::Arc *arc) {
    if(arc->event == nullptr)
        return;
    auto event = arc->event;
    if(event->next != nullptr)
        event->next->prev = event->prev;
    if(event->prev != nullptr)
        event->prev->next = event->next;
    else
        head = event->next;
    delete event;
    arc->event = nullptr;
}

bool FortuneAlgorithm::EventQueue::empty() {
    return head == nullptr;
}

FortuneAlgorithm::Event *FortuneAlgorithm::EventQueue::pop() {
    auto event = head;
    head = head->next;
    if(head != nullptr)
        head->prev = nullptr;
    else
        tail = nullptr;
    return event;
}

FortuneAlgorithm::EventQueue::~EventQueue() {
    while(head != nullptr) {
        auto next = head->next;
        delete head;
        head = next;
    }
}
