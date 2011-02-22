/**************************************************************************
 * FILE: red-black-tree.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 08-September-2009 
 * PROJECT: shared
 * COPYRIGHT: UQ, 2009
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: A red-black semi-balanced binary search tree. Usable as a set
 * or a map. Provides O(log(N)) lookup, insert and delete operations.
 **************************************************************************/
#include <assert.h>
#include <strings.h>
#include "red-black-tree.h"

/*
 * Structure for a red-black semi-balanced binary search tree.
 * A red-black tree has the following properties (from Wikipedia)
 * 1. A node is either red or black.
 * 2. The root is black.
 * 3. All leaves are black (in this case leaves are the null pointers)
 * 4. Both children of every red node are black.
 * 5. Every simple path from a given node to any of its descendant leaves
 *    contains the same number of black nodes.
 *
 * The result of these properties is that longest route to a leaf is at most
 * twice the shortest route. Lookup, Insert and Delete are all O(log(n)).
 */
struct rbtree_t {
  RBNODE_T *root;                     // the root node of the tree
  int size;                           // the number of internal nodes in the tree
  int (*key_compare) (void*, void*);  // a function to compare two keys
  void* (*key_copy)(void*);           // OPTIONAL a function to copy a key
  void (*key_free) (void*);           // OPTIONAL a function to free a key
  void* (*value_copy) (void*);        // OPTIONAL a function to copy a value
  void (*value_free) (void*);         // OPTIONAL a function to free a value
};

/*
 * Internal node in the red-black tree. Leaves are defined as existing vitually 
 * off any NULL left or right pointers. A NULL parent pointer implys that the node
 * is the root of the tree.
 */
struct rbnode_t {
  BOOLEAN_T is_red; // true if the node is a red node
  RBNODE_T *left;   // a child tree that has nodes with smaller keys than this one
  RBNODE_T *right;  // a child tree that has nodes with larger keys than this one
  RBNODE_T *parent; // the parent of this node in the tree
  void *key;        // the key of this node
  void *value;      // the value of this node
};

/*
 * check_recursive
 * Debugging function. Checks that the node and all subnodes are consistant with a red-black tree's structure.
 * Should only be called by rbtree_check.
 */
static int check_recursive(RBNODE_T *node, BOOLEAN_T must_be_black, int (*cmp)(void*, void*), int *black_nodes_to_leaf) {
  if (node->is_red && must_be_black) die("A node that must be black is red\n");
  int count = 1, left_black_nodes = 0, right_black_nodes = 0;
  if (node->left != NULL) {
    if (node->left->parent != node) die("Left node has wrong parent node\n"); 
    if (cmp(node->key, node->left->key) < 0) die("Left node has larger key\n");
    count += check_recursive(node->left, node->is_red, cmp, &left_black_nodes);
  }
  if (node->right != NULL) {
    if (node->right->parent != node) die("Right node has wrong parent node\n"); 
    if (cmp(node->key, node->right->key) > 0) die("Right node has smaller key\n");
    count += check_recursive(node->right, node->is_red, cmp, &right_black_nodes);
  }
  if (left_black_nodes != right_black_nodes)
    die("Number of black nodes in a simple path to a left leaf node must be the same as the right leaf node\n");
  *black_nodes_to_leaf = left_black_nodes;
  if (!(node->is_red)) *black_nodes_to_leaf += 1;
  return count;
}

/*
 * rbtree_check
 * Debugging function. Checks that the tree is consistant with a red-black tree's structure.
 */
void rbtree_check(RBTREE_T *tree) {
  if (tree == NULL) die("Tree is null\n");
  if (tree->key_compare == NULL) die("key_compare is null\n");
  if (tree->size == 0) {
    if (tree->root != NULL) die("Root expected to be null as tree is empty\n");
  } else {
    if (tree->root->parent != NULL) die("Root node has parent\n");
    int black_nodes_to_leaf = 0;
    int count = check_recursive(tree->root, TRUE, tree->key_compare, &black_nodes_to_leaf);
    if (count != tree->size)
      die("Mismatch between recorded size and actual node count\n");
  }
}

/*
 * rbtree_create
 * create a red-black tree.
 * takes a key comparator and optional key/value copy/free functions
 */
RBTREE_T *rbtree_create(
    int (*key_compare) (void*, void*), 
    void* (*key_copy)(void*), 
    void (*key_free) (void*), 
    void* (*value_copy) (void*), 
    void (*value_free) (void*)
  ) 
{
  RBTREE_T *tree = (RBTREE_T*)mm_malloc(sizeof(RBTREE_T));
  tree->root = NULL;
  tree->size = 0;
  tree->key_compare = key_compare;
  tree->key_copy = key_copy;
  tree->key_free = key_free;
  tree->value_copy = value_copy;
  tree->value_free = value_free;
  return tree;
}

/*
 * rbtree_destroy
 * destroys a red-black tree.
 * If any key or value free functions have been set then it will run 
 * them on the keys and values.
 */
void rbtree_destroy(RBTREE_T *tree) {
  RBNODE_T *current, *destroy;
  int destroyed = 0;
  current = tree->root;
  while (current != NULL) {
    // try changing to the left node, if that's not avaliable try the right
    // and if that's not avaliable then destroy the current node and change to 
    // the parent.
    if (current->left != NULL) {
      current = current->left;
      continue;
    }
    if (current->right != NULL) {
      current= current->right;
      continue;
    }
    destroy = current;
    current = current->parent;
    //remove references in the parent to the node being destroyed
    if (current != NULL) {
      if (current->left == destroy) {
        current->left = NULL;
      } else {
        current->right = NULL;
      }
    }
    // destroy the content of the node
    if (tree->key_free) tree->key_free(destroy->key);
    if (tree->value_free) tree->value_free(destroy->value);
    // destroy the node
    memset(destroy, 0, sizeof(RBNODE_T));
    free(destroy);
    ++destroyed;
  }
  assert(destroyed == tree->size);
  // destroy the tree
  memset(tree, 0, sizeof(RBTREE_T));
  free(tree);
}

/*
 * rbtree_alter_key_copy
 * Changes the function used to copy the key. Might be useful if you are
 * loading from multiple sources one of which you want to copy the key
 * and the other you just want to use it as is. If NULL is passed then
 * no attempt will be made to copy the keys.
 */
void rbtree_alter_key_copy(RBTREE_T *tree, void* (*key_copy)(void*)) {
  tree->key_copy = key_copy;
}

/*
 * rbtree_alter_key_free
 * Changes the function used to free the key.
 * If NULL is passed then no attempt will be made to free the keys.
 */
void rbtree_alter_key_free(RBTREE_T *tree, void (*key_free)(void*)) {
  tree->key_free = key_free;
}

/*
 * rbtree_alter_value_copy
 * Changes the function used to copy the value. Might be useful if you are
 * loading from multiple sources, one of which you want to copy the value
 * and the other you just want to use it as is. If NULL is passed then
 * no attempt will be made to copy the values.
 */
void rbtree_alter_value_copy(RBTREE_T *tree, void* (*value_copy)(void*)) {
  tree->value_copy = value_copy;
}

/*
 * rbtree_alter_value_free
 * Changes the function used to free the value. If NULL is passed then
 * no attempt will be made to free the values.
 */
void rbtree_alter_value_free(RBTREE_T *tree, void (*value_free)(void*)) {
  tree->value_free = value_free;
}

/*
 * rbtree_size
 * number of nodes in a red-black tree
 */
int rbtree_size(RBTREE_T *tree) {
  return tree->size;
}


/*
 * rotate_tree_right - rotates the tree right around the passed node
 * the passed node is moved down and to the right, while its left
 * branch is moved up and to the right.
 *      x                         y
 *     / \       right rotate    / \
 *    y   C      ---------->    A   x
 *   / \         (around x)        / \
 *  A   B                         B   C
 */
static inline void rotate_tree_right(RBTREE_T *tree, RBNODE_T *x) {
  RBNODE_T *y;
  assert(x != NULL);
  y = x->left;
  assert(y != NULL);
  y->parent = x->parent;
  if (x->parent) {
    if (x->parent->left == x) {
      x->parent->left = y;
    } else {
      x->parent->right = y;
    }
  } else {
    tree->root = y;
  }
  x->parent = y;
  x->left = y->right;
  if (x->left) x->left->parent = x;
  y->right = x;
}

/*
 * rotate_tree_left - rotates the tree left around the passed node.
 * the passed node is moved down and to the left, while its right
 * branch is moved up and to the left.
 *      x                         y
 *     / \       left rotate     / \
 *    A   y      ---------->    x   C
 *       / \     (around x)    / \
 *      B   C                 A   B
 */
static inline void rotate_tree_left(RBTREE_T *tree, RBNODE_T *x) {
  RBNODE_T *y;
  assert(x != NULL);
  y = x->right;
  assert(y != NULL);
  y->parent = x->parent;
  if (x->parent) {
    if (x->parent->left == x) {
      x->parent->left = y;
    } else {
      x->parent->right = y;
    }
  } else {
    tree->root = y;
  }
  x->parent = y;
  x->right = y->left;
  if (x->right) x->right->parent = x;
  y->left = x;
}


/*
 * fix_consecutive_reds
 * fixes consecutive reds in the red-black tree structure after an insert has occured.
 */
static inline void fix_consecutive_reds(RBTREE_T *tree, RBNODE_T *node) {
  RBNODE_T *current, *parent, *uncle, *grandparent, *greatgrandparent;
  current = node;
  parent = node->parent; // assumes this is non-null
  while (TRUE) {
    grandparent = parent->parent;
    if (grandparent == NULL) {
      // parent must be root node, we have fixed everything below so exit
      parent->is_red = FALSE;
      return;
    }
    // find the opposing branch (uncle)
    if (grandparent->left == parent) uncle = grandparent->right;
    else uncle = grandparent->left;

    if (uncle != NULL && uncle->is_red) {
      parent->is_red = FALSE;
      uncle->is_red = FALSE;
      greatgrandparent = grandparent->parent;
      if (greatgrandparent != NULL) {
        grandparent->is_red = TRUE;
        if (greatgrandparent->is_red) {
          current = grandparent;
          parent = greatgrandparent;
          continue;
        }
      }
    } else if (current == parent->left) {
      if (parent == grandparent->left) {
        rotate_tree_right(tree, grandparent);
        parent->is_red = FALSE;
        grandparent->is_red = TRUE;
      } else {
        rotate_tree_right(tree, parent);
        rotate_tree_left(tree, grandparent);
        current->is_red = FALSE;
        grandparent->is_red = TRUE;
      }
    } else {
      if (parent == grandparent->left) {
        rotate_tree_left(tree, parent);
        rotate_tree_right(tree, grandparent);
        current->is_red = FALSE;
        grandparent->is_red = TRUE;
      } else {
        rotate_tree_left(tree, grandparent);
        parent->is_red = FALSE;
        grandparent->is_red = TRUE;
      }
    }
    return;
  }
}

/*
 * rbtree_lookup
 * lookup a node in the tree. If the node doesn't exist and create is true then a new
 * node will be created using the key. If created is non-null then it will be set to true
 * when a new node is created.
 */
RBNODE_T* rbtree_lookup(RBTREE_T *tree, void *key, BOOLEAN_T create, BOOLEAN_T *created) {
  RBNODE_T *current, *next;
  BOOLEAN_T left = FALSE;
  int cmp;
  current = NULL;
  left = FALSE; //stop a compilier warning (though left is always initilized before use anyway)
  next = tree->root;
  while (next != NULL) {
    current = next;
    cmp = tree->key_compare(key, current->key);
    if (cmp == 0) { // found the key!
      if (created != NULL) *created = FALSE;
      return current;
    } else if (cmp < 0) {
      left = TRUE;
      next = current->left;
    } else {
      left = FALSE;
      next = current->right;
    }
  }
  if (create) {
    // create a new node
    next = mm_malloc(sizeof(RBNODE_T));
    next->left = NULL;
    next->right = NULL;
    next->value = NULL;
    if (tree->key_copy != NULL) {
      next->key = tree->key_copy(key);
    } else {
      next->key = key;
    }
    if (current == NULL) { // new root
      next->is_red = FALSE;
      next->parent = NULL;
      tree->root = next;
    } else { // new internal node
      next->is_red = TRUE;
      next->parent = current;
      if (left) current->left = next;
      else current->right = next;
      if (current->is_red) fix_consecutive_reds(tree, next);
    }
    tree->size += 1;
    // set the created output value
    if (created != NULL) *created = TRUE;
    return next;
  }
  if (created != NULL) *created = FALSE;
  return NULL;
}

/*
 * rbtree_get
 * Gets the value for a key. If the key doesn't exist then NULL is returned.
 */
void *rbtree_get(RBTREE_T *tree, void *key) {
  RBNODE_T *node;
  node = rbtree_lookup(tree, key, FALSE, NULL);
  if (node == NULL) return NULL;
  return node->value;
}

/*
 * rbnode_get
 * Gets the value from a node.
 */
void *rbnode_get(RBNODE_T *node) {
	if (node == NULL) return NULL;
	return node->value;
}

/*
 * rbtree_set
 * Updates the value. If the new value equals the previous value by pointer comparison then
 * nothing is done. Otherwise, if there was a previous value it will be freed using value_free 
 * (if it was passed to the constructor) and the new value will be copied using value_copy
 * (again only if it was passed to the constructor).
 */
void rbtree_set(RBTREE_T *tree, RBNODE_T *node, void *value) {
  if (node->value != value) {
    if (node->value) {
      if (tree->value_free) tree->value_free(node->value); 
    }
    if (tree->value_copy) {
      node->value = tree->value_copy(value);
    } else {
      node->value = value;
    }
  }
}

/*
 * rbtree_put
 * Puts a key, value combination into the tree. If the key already exists the value is set to
 * the passed value as described in rbtree_set. Returns the node in the tree that represents
 * the key, value combination.
 */
RBNODE_T* rbtree_put(RBTREE_T *tree, void *key, void *value) {
  RBNODE_T *node;
  node = rbtree_lookup(tree, key, TRUE, NULL);
  rbtree_set(tree, node, value);
  return node;
}

/*
 * fix_double_black - rebalances a red-black tree before a delete so that the 
 * deletion won't effect the black height of a tree.
 * assumes that double_black_node is a double black node, which is a virtual
 * state a black node can be in when it combines with it's leaf node just
 * prior to deletion.
 */
static inline void fix_double_black(RBTREE_T *tree, RBNODE_T *double_black_node) {
  RBNODE_T *current, *parent, *sibling, *far_nephew, *near_nephew;
  current = double_black_node;
  parent = double_black_node->parent;

  while (TRUE) {
    if (current == parent->left) { // current is the left child of parent
      sibling = parent->right;
      if (sibling->is_red) {
        // the sibling is red, which we don't want... rearrange the structure so it's black
        parent->is_red = TRUE;
        sibling->is_red = FALSE;
        rotate_tree_left(tree, parent);
        sibling = parent->right;
      }
      far_nephew = sibling->right;
      near_nephew = sibling->left;
      if (far_nephew && far_nephew->is_red) {
        rotate_tree_left(tree, parent);
        sibling->is_red = parent->is_red;
        parent->is_red = FALSE;
        far_nephew->is_red = FALSE;
        return;
      } else if (near_nephew && near_nephew->is_red) {
        rotate_tree_right(tree, sibling);
        rotate_tree_left(tree, parent);
        near_nephew->is_red = parent->is_red;
        parent->is_red = FALSE;
        return;
      }
    } else { // current is the right child of parent
      sibling = parent->left;
      if (sibling->is_red) {
        // the sibling is red, which we don't want... rearrange the structure so it's black
        parent->is_red = TRUE;
        sibling->is_red = FALSE;
        rotate_tree_right(tree, parent);
        sibling = parent->left;
      }
      far_nephew = sibling->left;
      near_nephew = sibling->right;
      if (far_nephew && far_nephew->is_red) {
        rotate_tree_right(tree, parent);
        sibling->is_red = parent->is_red;
        parent->is_red = FALSE;
        far_nephew->is_red = FALSE;
        return;
      } else if (near_nephew && near_nephew->is_red) {
        rotate_tree_left(tree, sibling);
        rotate_tree_right(tree, parent);
        near_nephew->is_red = parent->is_red;
        parent->is_red = FALSE;
        return;
      }
    }
    // neither nephew is red so both must be black
    sibling->is_red = TRUE;
    if (parent->is_red) {
      parent->is_red = FALSE;
      return;
    }
    //parent is now double black as it was black already!
    current = parent;
    parent = current->parent;
    if (parent == NULL) {
      // found the root node
      // as double black means nothing for the root then just return
      return;
    }
  }
}

/*
 * delete_internal
 * Removes a passed node from the passed tree and ensures the tree maintains
 * the red-black tree properties. Does not alter the node.
 */
static void delete_internal(RBTREE_T *tree, RBNODE_T *node) {
  RBNODE_T *left, *right, *parent, *successor;

  parent = node->parent;
  left = node->left;
  right = node->right;

  if (left == NULL && right == NULL) {
    if (parent == NULL) { // deleting root node
      tree->root = NULL;
      tree->size -= 1;
      assert(tree->size == 0);
      return;
    }
    if (!(node->is_red)) { // node is black, deleting causes double black!
      fix_double_black(tree, node);
    }
    // remove reference to the node from the tree
    if (parent->left == node) { 
      parent->left = NULL;
    } else {
      parent->right = NULL;
    }
  } else if (right == NULL) { // left subnode exists, move it up to replace
    if (parent == NULL) tree->root = left;
    else if (parent->left == node) parent->left = left;
    else parent->right = left;
    left->parent = parent;
    if (!(node->is_red)) left->is_red = FALSE;
  } else if (left == NULL) { // right subnode exists, move it up to replace
    if (parent == NULL) tree->root = right;
    else if (parent->left == node) parent->left = right;
    else parent->right = right;
    right->parent = parent;
    if (!(node->is_red)) right->is_red = FALSE;
  } else {
    //note that when we only had one child node we could be certain
    //that the balance wasn't made worse by simply swapping it in but
    //with both, we have to delete the successor node and use it as a 
    //replacement
    
    //get the next larger node (successor) 
    successor = right;
    while (successor->left) successor = successor->left;

    // delete the successor from the tree where it is currently
    delete_internal(tree, successor);
    // update right, left and parent
    right = node->right;
    left = node->left;
    parent = node->parent;
    // replace node with successor
    successor->left = left;
    successor->right = right;
    successor->parent = parent;
    if (parent == NULL) tree->root = successor;
    else if (parent->left == node) parent->left = successor;
    else parent->right = successor;
    if (left) left->parent = successor;
    if (right) right->parent = successor;
    successor->is_red = node->is_red;
    return; // we don't want to deincrement size twice so exit first
  }
  tree->size -= 1;
}

/*
 * rbtree_delete
 * removes the passed node from the passed tree and optionally returns the removed content. 
 * If removed_key is non-null then the free_key (from the constructor) will not be run on the
 * key, but it will instead be returned in the removed_key pointer. If removed_value is 
 * non-null then free_value will not be run on the value, but it will instad be returned in
 * the removed_value pointer. The node will be destroyed.
 */
void rbtree_delete(RBTREE_T *tree, RBNODE_T *node, void **removed_key, void **removed_value) {
  delete_internal(tree, node);
  if (removed_key) {
    *removed_key = node->key;
  } else {
    if (tree->key_free) tree->key_free(node->key);
  }
  if (removed_value) {
    *removed_value = node->value;
  } else {
    if (tree->value_free) tree->value_free(node->value);
  }
  memset(node, 0, sizeof(RBNODE_T));
  free(node);
}

/*
 * rbtree_remove
 * Removes the key (and any value it has) from the tree. The key and value stored by the tree
 * are freed by free_key and free_value if these functions were set in the constructor.
 * Returns true if the key existed.
 */
BOOLEAN_T rbtree_remove(RBTREE_T *tree, void *key) {
  RBNODE_T *node;
  node = rbtree_lookup(tree, key, FALSE, NULL);
  if (node) {
    rbtree_delete(tree, node, NULL, NULL);
    return TRUE;
  }
  return FALSE;
}

/*
 * rbtree_first
 * Returns the smallest node in the tree
 */
RBNODE_T *rbtree_first(RBTREE_T *tree) {
  RBNODE_T *node;
  node = tree->root;
  if (node == NULL) return NULL;
  while (node->left) node = node->left;
  return node;
}

/*
 * rbtree_last
 * Returns the largest node in the tree
 */
RBNODE_T *rbtree_last(RBTREE_T *tree) {
  RBNODE_T *node;
  node = tree->root;
  if (node == NULL) return NULL;
  while (node->right) node = node->right;
  return node;
}

/*
 * rbtree_next
 * Returns the next larger node in the tree or null if node is the largest
 */
RBNODE_T *rbtree_next(RBNODE_T *node) {
  RBNODE_T *parent;
  if (node->right == NULL) {
    //no right branch so find the first parent that has this tree as it's left branch
    parent = node->parent;
    while (parent != NULL && parent->right == node) {
      node = parent;
      parent = node->parent;
    }
    // node is the left tree of parent so parent is next larger
    return parent;
  } else {
    //has right branch, so next is the left-most on the right
    node = node->right;
    while (node->left) node = node->left;
    return node;
  }
}

/*
 * rbtree_prev
 * Returns the previous smaller node in the tree or null if node is the smallest
 */
RBNODE_T *rbtree_prev(RBNODE_T *node) {
  RBNODE_T *parent;
  if (node->left == NULL) {
    parent = node->parent;
    while (parent != NULL && parent->left == node) {
      node = parent;
      parent = node->parent;
    }
    // node is the right tree of parent so parent is next smaller
    return parent;
  } else {
    //hash left branch, so previous is right-most on the left
    node = node->left;
    while (node->right) node = node->right;
    return node;
  }
}

/*
 * rbtree_key
 * Returns the key. This key is used internally by the red black tree so do not
 * modify it as it will cause undefined results.
 */
void *rbtree_key(RBNODE_T *node) {
  return node->key;
}

/*
 * rbtree_value
 * Returns the value.
 */
void *rbtree_value(RBNODE_T *node) {
  return node->value;
}


/*
 * rbtree_strcmp
 * Utility function for using the red-black tree with normal strings.
 * Returns the result of comparing two strings using the strcmp function.
 */
int rbtree_strcmp(void *p1, void *p2) {
  return strcmp((char*)p1, (char*)p2);
}

/*
 * rbtree_strcpy
 * Utility function for using the red-black tree with normal strings.
 * Note that the free function can be used as the counterpart.
 * Returns a copy of the passed string.
 */
void* rbtree_strcpy(void *p) {
  char *original = (char*)p;
  int size = (strlen(original) + 1);
  char *copy = (char*)mm_malloc(sizeof(char) * size);
  strncpy(copy, original, size);
  return copy;
}

/*
 * rbtree_intcmp
 * Utility function for using the red-black tree with pointers to integers.
 * Returns -1, 0 or 1 respectively if ints p1 < p2, p1 == p2 or p1 > p2.
 */
int rbtree_intcmp(void *p1, void *p2) {
  int i1, i2;
  i1 = *((int*)p1);
  i2 = *((int*)p2);
  if (i1 < i2) {
    return -1;
  } else if (i1 == i2) {
    return 0;
  } else {
    return 1;
  }
}

/*
 * rbtree_intcpy
 * Utility function for using the red-black tree with pointers to integers.
 * Note that the free function can be used as the counterpart.
 * Returns a malloc'ed copy of the passed int.
 */
void* rbtree_intcpy(void *p) {
  int *copy = (int*)mm_malloc(sizeof(int));
  *copy = *((int*)p);
  return copy;
}
