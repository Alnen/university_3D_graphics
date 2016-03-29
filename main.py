from vispy import app, visuals, scene, geometry, util, gloo
from vispy.app.canvas import DrawEvent
from vispy.gloo.util import _screenshot
from moviepy.editor import VideoClip
import numpy as np
import math
import typing


class Canvas(scene.SceneCanvas):
    def __init__(self, transform_function, *args, **kwargs):
        scene.SceneCanvas.__init__(self, *args, **kwargs)
        self.timer = app.Timer('auto', connect=self.update, start=True)
        self.view = self.central_widget.add_view()
        self.view.camera = scene.cameras.TurntableCamera()
        self.transform_function = transform_function
        self.show()

    def on_draw(self, event):
        self.transform_function(self.view.scene.children)
        scene.SceneCanvas.on_draw(self, event)

    def animation(self, t):
        self.on_draw(DrawEvent('draw'))
        return _screenshot((0, 0, self.size[0], self.size[1]))[:, :, :3]


def get_rotation_matrix(*, start_point: np.ndarray, end_point: typing.Any) -> np.ndarray:
    if type(end_point) is tuple:
        end_point = np.array(end_point, np.float32)
    if np.array_equal(start_point, end_point):
        return np.identity(3, np.float32)
    if np.array_equal(start_point, -end_point):
        if start_point[0] == 1.0:
            return util.transforms.rotate(180, (0, 0, 1))[:3, :3]
        else:
            return util.transforms.rotate(180, (1, 0, 0))[:3, :3]
    else:
        x = np.cross(start_point, end_point)
        x = np.divide(x, np.linalg.norm(x))
        start_point_length = np.linalg.norm(start_point)
        end_point_length = np.linalg.norm(end_point)
        cos_theta = np.dot(start_point, end_point)/(start_point_length*end_point_length)
        theta = math.acos(cos_theta)

        A = np.array([[    0, -x[2],  x[1]],
                      [ x[2],     0, -x[0]],
                      [-x[1],  x[0],    0]],
                     np.float32)
        rotate_matrix = np.identity(3, np.float32)
        rotate_matrix += A*math.sin(theta)
        rotate_matrix += np.dot(A, A)*(1 - math.cos(theta))
        if x[2] <= 0:
            rotate_matrix = np.dot(rotate_matrix, util.transforms.rotate(180, (0, 0, 1))[:3, :3])
        return rotate_matrix


def generate_circle(segments, radius=1, normal=(0, 0, 1)):
    angle = (2 * math.pi) / segments
    vertex_buffer = np.array([(radius*math.cos(i*angle), radius*math.sin(i*angle), 0)
                              for i in range(segments)], np.float32)
    vertex_buffer = np.dot(vertex_buffer, get_rotation_matrix(start_point=(0, 0, 1), end_point=normal))
    return vertex_buffer


def truncated_icosahedron_data():
    sphere_radius = 0.4
    sphere_rows = 10
    sphere_cols = 10
    cylinder_radius = 0.25 * sphere_radius
    cylinder_rows = 10
    cylinder_cols = 10
    color = (1.0, 0.725, 0.275)

    truncated_icosahedron_vertices_coordinates = np.array(
        [
            (0, 1, 4.854101966249685),
            (0, 1, -4.854101966249685),
            (0, -1, 4.854101966249685),
            (0, -1, -4.854101966249685),
            (1, 4.854101966249685, 0),
            (1, -4.854101966249685, 0),
            (-1, 4.854101966249685, 0),
            (-1, -4.854101966249685, 0),
            (4.854101966249685, 0, 1),
            (4.854101966249685, 0, -1),
            (-4.854101966249685, 0, 1),
            (-4.854101966249685, 0, -1),
            (2, 4.23606797749979, 1.618033988749895),
            (2, 4.23606797749979, -1.618033988749895),
            (2, -4.23606797749979, 1.618033988749895),
            (2, -4.23606797749979, -1.618033988749895),
            (-2, 4.23606797749979, 1.618033988749895),
            (-2, 4.23606797749979, -1.618033988749895),
            (-2, -4.23606797749979, 1.618033988749895),
            (-2, -4.23606797749979, -1.618033988749895),
            (4.23606797749979, 1.618033988749895, 2),
            (4.23606797749979, 1.618033988749895, -2),
            (4.23606797749979, -1.618033988749895, 2),
            (4.23606797749979, -1.618033988749895, -2),
            (-4.23606797749979, 1.618033988749895, 2),
            (-4.23606797749979, 1.618033988749895, -2),
            (-4.23606797749979, -1.618033988749895, 2),
            (-4.23606797749979, -1.618033988749895, -2),
            (1.618033988749895, 2, 4.23606797749979),
            (1.618033988749895, 2, -4.23606797749979),
            (1.618033988749895, -2, 4.23606797749979),
            (1.618033988749895, -2, -4.23606797749979),
            (-1.618033988749895, 2, 4.23606797749979),
            (-1.618033988749895, 2, -4.23606797749979),
            (-1.618033988749895, -2, 4.23606797749979),
            (-1.618033988749895, -2, -4.23606797749979),
            (1, 3.618033988749895, 3.23606797749979),
            (1, 3.618033988749895, -3.23606797749979),
            (1, -3.618033988749895, 3.23606797749979),
            (1, -3.618033988749895, -3.23606797749979),
            (-1, 3.618033988749895, 3.23606797749979),
            (-1, 3.618033988749895, -3.23606797749979),
            (-1, -3.618033988749895, 3.23606797749979),
            (-1, -3.618033988749895, -3.23606797749979),
            (3.618033988749895, 3.23606797749979, 1),
            (3.618033988749895, 3.23606797749979, -1),
            (3.618033988749895, -3.23606797749979, 1),
            (3.618033988749895, -3.23606797749979, -1),
            (-3.618033988749895, 3.23606797749979, 1),
            (-3.618033988749895, 3.23606797749979, -1),
            (-3.618033988749895, -3.23606797749979, 1),
            (-3.618033988749895, -3.23606797749979, -1),
            (3.23606797749979, 1, 3.618033988749895),
            (3.23606797749979, 1, -3.618033988749895),
            (3.23606797749979, -1, 3.618033988749895),
            (3.23606797749979, -1, -3.618033988749895),
            (-3.23606797749979, 1, 3.618033988749895),
            (-3.23606797749979, 1, -3.618033988749895),
            (-3.23606797749979, -1, 3.618033988749895),
            (-3.23606797749979, -1, -3.618033988749895)
        ],
        np.float32)

    edges_between_vertices = [
        (0, 2), (4, 6), (8, 9), (12, 4), (0, 32), (32, 56), (56, 58), (34, 58), (34, 2), (34, 42), (42, 18),
        (18, 50), (50, 26), (26, 58), (18, 7), (7, 19), (19, 51), (50, 51), (7, 5), (5, 15), (15, 39), (39, 43),
        (43, 19), (15, 47), (47, 23), (55, 31), (31, 39), (23, 55), (23, 9), (9, 21), (8, 22), (22, 46), (46, 47),
        (5, 14), (14, 46), (8, 20), (20, 44), (44, 45), (45, 21), (55, 53), (53, 21), (44, 12), (4, 13), (45, 13),
        (13, 37), (37, 29), (29, 53), (14, 38), (38, 42), (38, 30), (30, 2), (30, 54), (54, 22), (54, 52), (52, 20),
        (0, 28), (28, 52), (28, 36), (36, 12), (36, 40), (40, 16), (16, 6), (6, 17), (17, 41), (41, 37), (56, 24),
        (24, 10), (10, 26), (10, 11), (11, 27), (27, 51), (40, 32), (16, 48), (48, 24), (48, 49), (49, 25), (25, 11),
        (25, 57), (57, 59), (59, 27), (59, 35), (35, 43), (17, 49), (41, 33), (33, 1), (1, 29), (35, 3), (3, 31),
        (1, 3), (33, 57)
    ]

    return (truncated_icosahedron_vertices_coordinates,
            edges_between_vertices,
            sphere_radius,
            cylinder_radius,
            sphere_rows,
            sphere_cols,
            cylinder_rows,
            cylinder_cols,
            color)


def torus_data(rows=30, cols=15):
    sphere_radius = 0.1
    sphere_rows = 10
    sphere_cols = 10
    cylinder_radius = 0.25 * sphere_radius
    cylinder_rows = 10
    cylinder_cols = 10
    color = (0.1, 0.6, 0.1)

    def torus_vertices_coordinates_generator():
        torus_global_radius = 8
        angle_local = 180 / cols
        angle_global = 360 / rows
        circle_vertices_coordinates = generate_circle(cols, 2, (1, 0, 0))
        shifted_circle_vertices_coordinates = np.dot(
            circle_vertices_coordinates,
            util.transforms.rotate(angle_local, (1, 0, 0))[:3, :3])
        circle_vertices_coordinates += (0, torus_global_radius, 0)
        shifted_circle_vertices_coordinates += (0, torus_global_radius, 0)
        shifted_circle_vertices_coordinates = np.dot(
            shifted_circle_vertices_coordinates,
            util.transforms.rotate(angle_global/3, (0, 0, 1))[:3, :3])

        for row in range(rows):
            if row % 2 == 0:
                rotation_matrix = util.transforms.rotate(row*angle_global, (0, 0, 1))[:3, :3]
                for vertex in np.dot(circle_vertices_coordinates, rotation_matrix):
                    yield vertex
                for vertex in np.dot(shifted_circle_vertices_coordinates, rotation_matrix):
                    yield vertex
            else:
                rotation_matrix = util.transforms.rotate(row*angle_global - angle_global/3, (0, 0, 1))[:3, :3]
                for vertex in np.dot(shifted_circle_vertices_coordinates, rotation_matrix):
                    yield vertex

                rotation_matrix = util.transforms.rotate(row*angle_global + angle_global/3, (0, 0, 1))[:3, :3]
                for vertex in np.dot(circle_vertices_coordinates, rotation_matrix):
                    yield vertex

    def edges_between_vertices():
        for i in range(0, 2*rows, 2):
            for j in range(cols):
                first_idx = i * cols + j
                second_idx = ((i + 1) * cols) + j
                if i % 4 == 0:
                    third_idx = ((i + 1) * cols) + (j - 1) % cols
                else:
                    third_idx = ((i + 1) * cols) + (j + 1) % cols
                fourth_idx = ((i - 1) * cols) % (2 * rows * cols) + j
                yield first_idx, second_idx
                yield first_idx, third_idx
                yield first_idx, fourth_idx

    return (list(torus_vertices_coordinates_generator()),
            edges_between_vertices(),
            sphere_radius,
            cylinder_radius,
            sphere_rows,
            sphere_cols,
            cylinder_rows,
            cylinder_cols,
            color)


def generate_skeleton(
        truncated_icosahedron_vertices_coordinates,
        edges_between_vertices,
        sphere_radius=0.4,
        cylinder_radius=0.1,
        sphere_rows=10,
        sphere_cols=10,
        cylinder_rows=10,
        cylinder_cols=10,
        color=(1, 0.6, 0)):
    spheres_geometry = []
    for i, vertex_coordinate in enumerate(truncated_icosahedron_vertices_coordinates):
        new_mesh = geometry.create_sphere(sphere_rows, sphere_cols, radius=sphere_radius)
        new_mesh.set_vertices(new_mesh.get_vertices() + vertex_coordinate, reset_normals=False)
        spheres_geometry.append(new_mesh)

    cylinders_geometry = []
    for first_vertex_index, second_vertex_index in edges_between_vertices:
        first_vertex_coordinate = truncated_icosahedron_vertices_coordinates[first_vertex_index]
        second_vertex_coordinate = truncated_icosahedron_vertices_coordinates[second_vertex_index]
        cylinder = generate_cylinder_geometry(first_vertex_coordinate, second_vertex_coordinate,
                                              cylinder_rows, cylinder_cols, cylinder_radius, sphere_radius)
        cylinders_geometry.append(cylinder)

    final_object_vertices, final_object_faces = append_geometries(np.empty([0, 3], np.float32),
                                                                  np.empty([0, 3], np.uint32),
                                                                  spheres_geometry)
    final_object_vertices, final_object_faces = append_geometries(final_object_vertices, final_object_faces,
                                                                  cylinders_geometry)

    return scene.visuals.Mesh(vertices=final_object_vertices, faces=final_object_faces, color=color, shading='smooth')


def append_geometries(vertices, faces, geometries):
    for geometry in geometries:
        faces = np.append(faces, geometry.get_faces() + vertices.shape[0], axis=0)
        vertices = np.append(vertices, geometry.get_vertices(), axis=0)
    return vertices, faces


def generate_cylinder_geometry(first_vertex_coordinate, second_vertex_coordinate,
                               cylinder_rows, cylinder_cols, cylinder_radius, sphere_radius):
    middle_point = (first_vertex_coordinate + second_vertex_coordinate) / 2
    distance_between_vertices = np.linalg.norm(first_vertex_coordinate - second_vertex_coordinate)
    normalized_vector = (second_vertex_coordinate - first_vertex_coordinate) / distance_between_vertices
    cylinder_length = distance_between_vertices - 1.75 * sphere_radius
    cylinder = geometry.create_cylinder(
        cylinder_rows, cylinder_cols, radius=[cylinder_radius, cylinder_radius], length=cylinder_length)
    vertices = cylinder.get_vertices()
    vertices += (0, 0, -cylinder_length / 2)
    rotation_matrix = get_rotation_matrix(start_point=np.array((0, 0, 1), np.float32), end_point=normalized_vector)
    vertices = np.dot(vertices, rotation_matrix)
    vertices += middle_point
    cylinder.set_vertices(vertices)
    return cylinder


def generate_move_rotate_transform_chain_wrapper(wrapped_object):
    wrapper = scene.Node()
    wrapper.transform = visuals.transforms.ChainTransform([
        visuals.transforms.AffineTransform(),
        visuals.transforms.AffineTransform()
    ])
    wrapped_object.parent = wrapper
    return wrapper


class Rotator:
    def __init__(self, object_normal_list, angle_in_second, expected_framerate=60):
        self.object_normal_list = object_normal_list
        self.velocity_per_frame = angle_in_second/expected_framerate

    def transform(self, objects):
        for index, normal in self.object_normal_list:
            objects[index].transform.transforms[1].rotate(self.velocity_per_frame, normal)


def generate_pov_ray_script(data_generators):
    def vertex_to_string(elem):
        return "    <{}, {}, {}>".format(elem[0], elem[1], elem[2])

    def edge_to_string(elem):
        return "    <{}, {}>".format(elem[0], elem[1])

    def transform_verticies(vertices, transform_matrix):
        for vertex in vertices:
            yield np.dot(vertex, transform_matrix)
    array_declaration_template = "#declare {} = array[{}]\n{{\n{}\n}}\n"
    variable_declaration_template = "#declare {} = {};\n"

    output = ""
    for data_generator, name, transform_matrix in data_generators:
        vertices, edges,  sphere_radius, cylinder_radius, _, _, _, _, _ = data_generator()
        output += variable_declaration_template.format("{}_sphere_radius".format(name), sphere_radius)
        output += variable_declaration_template.format("{}_cylinder_radius".format(name), cylinder_radius)

        vertices_generator = transform_verticies(vertices, transform_matrix)
        array_content, element_count = transform_data_to_povray_array_data(vertices_generator, vertex_to_string)
        output += variable_declaration_template.format("{}_vertices_count".format(name), element_count)
        output += array_declaration_template.format("{}_vertices".format(name), element_count, array_content)

        array_content, element_count = transform_data_to_povray_array_data(edges, edge_to_string)
        output += variable_declaration_template.format("{}_edges_count".format(name), element_count)
        output += array_declaration_template.format("{}_edges".format(name), element_count, array_content)

    return output


def transform_data_to_povray_array_data(vertices, element_to_string):
    count = 0

    def generator_element_counter(element):
        nonlocal count
        count += 1
        return element_to_string(element)

    array_content = ',\n'.join(generator_element_counter(vertex) for vertex in vertices)
    return array_content, count


def calculate_t_gold_2a_material():
    gold_base_color = np.array((1.00, 0.875, 0.575))
    diffuse_base_color = gold_base_color - np.array((0.00, 0.20, 0.40))
    color = gold_base_color - np.array((0.00, 0.15, 0.30))

    reflection_color = diffuse_base_color * 0.30 + np.array((0.25, 0.25, 0.25))
    ambient_color = diffuse_base_color * 0.12 + np.array((0.1, 0.1, 0.1))

    diffuse_coefficient = 1-(((reflection_color[0] + reflection_color[1] + reflection_color[2])/3)
                             + ((ambient_color[0] + ambient_color[1] + ambient_color[2])/3))
    diffuse_coefficient = max(diffuse_coefficient, 0)

    metallic = 1.0
    specular = 0.20
    roughness = 1.0/20.0

    return color, diffuse_coefficient, metallic, specular, roughness


if __name__ == '__main__':
    rotator = Rotator(object_normal_list=[(2, (1, 1, 0)), (3, (1, 0, 0))], angle_in_second=20)

    canvas = Canvas(rotator.transform, keys='interactive', show=True, app='pyqt5')

    truncated_icosahedron = generate_skeleton(*truncated_icosahedron_data())
    truncated_icosahedron = generate_move_rotate_transform_chain_wrapper(truncated_icosahedron)
    truncated_icosahedron.parent = canvas.view.scene

    torus = generate_skeleton(*torus_data(30, 15))
    torus = generate_move_rotate_transform_chain_wrapper(torus)
    torus.parent = canvas.view.scene

    pov_ray_arrays_declaration = generate_pov_ray_script([
        (truncated_icosahedron_data, 'truncated_icosahedron', util.transforms.rotate(70, (1, 0, 1))[:3, :3]),
        (torus_data, 'torus', util.transforms.rotate(70, (1, 0, 1))[:3, :3]@util.transforms.rotate(-70, (0, 1, 0))[:3, :3])
    ])

    output_file = open('lab2_data.inc', 'w')
    output_file.write(pov_ray_arrays_declaration)

    app.run()

    #clip = VideoClip(canvas.animation, duration=5).resize(0.3)
    #clip.write_gif('lab1.gif', fps=24, opt='OptimizePlus')
