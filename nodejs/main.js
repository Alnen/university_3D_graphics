var camera, scene, renderer;
var cameraControls;

function draw_scene() {
    init();
    animate();
}

function truncated_icosahedron_vertices_coordinates_data() {
    var sphere_radius = 40;
    var sphere_rows = 10;
    var sphere_cols = 10;
    var cylinder_radius = 0.25 * sphere_radius;
    var cylinder_rows = 10;
    var cylinder_cols = 10;
    var material = new THREE.MeshPhongMaterial();
    material.color = new THREE.Color(1.0, 0.725, 0.275);

    var truncated_icosahedron_vertices_coordinates = [
        new THREE.Vector3(0, 1, 4.854101966249685),
        new THREE.Vector3(0, 1, -4.854101966249685),
        new THREE.Vector3(0, -1, 4.854101966249685),
        new THREE.Vector3(0, -1, -4.854101966249685),
        new THREE.Vector3(1, 4.854101966249685, 0),
        new THREE.Vector3(1, -4.854101966249685, 0),
        new THREE.Vector3(-1, 4.854101966249685, 0),
        new THREE.Vector3(-1, -4.854101966249685, 0),
        new THREE.Vector3(4.854101966249685, 0, 1),
        new THREE.Vector3(4.854101966249685, 0, -1),
        new THREE.Vector3(-4.854101966249685, 0, 1),
        new THREE.Vector3(-4.854101966249685, 0, -1),
        new THREE.Vector3(2, 4.23606797749979, 1.618033988749895),
        new THREE.Vector3(2, 4.23606797749979, -1.618033988749895),
        new THREE.Vector3(2, -4.23606797749979, 1.618033988749895),
        new THREE.Vector3(2, -4.23606797749979, -1.618033988749895),
        new THREE.Vector3(-2, 4.23606797749979, 1.618033988749895),
        new THREE.Vector3(-2, 4.23606797749979, -1.618033988749895),
        new THREE.Vector3(-2, -4.23606797749979, 1.618033988749895),
        new THREE.Vector3(-2, -4.23606797749979, -1.618033988749895),
        new THREE.Vector3(4.23606797749979, 1.618033988749895, 2),
        new THREE.Vector3(4.23606797749979, 1.618033988749895, -2),
        new THREE.Vector3(4.23606797749979, -1.618033988749895, 2),
        new THREE.Vector3(4.23606797749979, -1.618033988749895, -2),
        new THREE.Vector3(-4.23606797749979, 1.618033988749895, 2),
        new THREE.Vector3(-4.23606797749979, 1.618033988749895, -2),
        new THREE.Vector3(-4.23606797749979, -1.618033988749895, 2),
        new THREE.Vector3(-4.23606797749979, -1.618033988749895, -2),
        new THREE.Vector3(1.618033988749895, 2, 4.23606797749979),
        new THREE.Vector3(1.618033988749895, 2, -4.23606797749979),
        new THREE.Vector3(1.618033988749895, -2, 4.23606797749979),
        new THREE.Vector3(1.618033988749895, -2, -4.23606797749979),
        new THREE.Vector3(-1.618033988749895, 2, 4.23606797749979),
        new THREE.Vector3(-1.618033988749895, 2, -4.23606797749979),
        new THREE.Vector3(-1.618033988749895, -2, 4.23606797749979),
        new THREE.Vector3(-1.618033988749895, -2, -4.23606797749979),
        new THREE.Vector3(1, 3.618033988749895, 3.23606797749979),
        new THREE.Vector3(1, 3.618033988749895, -3.23606797749979),
        new THREE.Vector3(1, -3.618033988749895, 3.23606797749979),
        new THREE.Vector3(1, -3.618033988749895, -3.23606797749979),
        new THREE.Vector3(-1, 3.618033988749895, 3.23606797749979),
        new THREE.Vector3(-1, 3.618033988749895, -3.23606797749979),
        new THREE.Vector3(-1, -3.618033988749895, 3.23606797749979),
        new THREE.Vector3(-1, -3.618033988749895, -3.23606797749979),
        new THREE.Vector3(3.618033988749895, 3.23606797749979, 1),
        new THREE.Vector3(3.618033988749895, 3.23606797749979, -1),
        new THREE.Vector3(3.618033988749895, -3.23606797749979, 1),
        new THREE.Vector3(3.618033988749895, -3.23606797749979, -1),
        new THREE.Vector3(-3.618033988749895, 3.23606797749979, 1),
        new THREE.Vector3(-3.618033988749895, 3.23606797749979, -1),
        new THREE.Vector3(-3.618033988749895, -3.23606797749979, 1),
        new THREE.Vector3(-3.618033988749895, -3.23606797749979, -1),
        new THREE.Vector3(3.23606797749979, 1, 3.618033988749895),
        new THREE.Vector3(3.23606797749979, 1, -3.618033988749895),
        new THREE.Vector3(3.23606797749979, -1, 3.618033988749895),
        new THREE.Vector3(3.23606797749979, -1, -3.618033988749895),
        new THREE.Vector3(-3.23606797749979, 1, 3.618033988749895),
        new THREE.Vector3(-3.23606797749979, 1, -3.618033988749895),
        new THREE.Vector3(-3.23606797749979, -1, 3.618033988749895),
        new THREE.Vector3(-3.23606797749979, -1, -3.618033988749895)
    ];

    var edges_between_vertices = [
        new THREE.Vector2(0, 2), new THREE.Vector2(4, 6), new THREE.Vector2(8, 9), new THREE.Vector2(12, 4), new THREE.Vector2(0, 32), new THREE.Vector2(32, 56), new THREE.Vector2(56, 58), new THREE.Vector2(34, 58), new THREE.Vector2(34, 2), new THREE.Vector2(34, 42), new THREE.Vector2(42, 18),
        new THREE.Vector2(18, 50), new THREE.Vector2(50, 26), new THREE.Vector2(26, 58), new THREE.Vector2(18, 7), new THREE.Vector2(7, 19), new THREE.Vector2(19, 51), new THREE.Vector2(50, 51), new THREE.Vector2(7, 5), new THREE.Vector2(5, 15), new THREE.Vector2(15, 39), new THREE.Vector2(39, 43),
        new THREE.Vector2(43, 19), new THREE.Vector2(15, 47), new THREE.Vector2(47, 23), new THREE.Vector2(55, 31), new THREE.Vector2(31, 39), new THREE.Vector2(23, 55), new THREE.Vector2(23, 9), new THREE.Vector2(9, 21), new THREE.Vector2(8, 22), new THREE.Vector2(22, 46), new THREE.Vector2(46, 47),
        new THREE.Vector2(5, 14), new THREE.Vector2(14, 46), new THREE.Vector2(8, 20), new THREE.Vector2(20, 44), new THREE.Vector2(44, 45), new THREE.Vector2(45, 21), new THREE.Vector2(55, 53), new THREE.Vector2(53, 21), new THREE.Vector2(44, 12), new THREE.Vector2(4, 13), new THREE.Vector2(45, 13),
        new THREE.Vector2(13, 37), new THREE.Vector2(37, 29), new THREE.Vector2(29, 53), new THREE.Vector2(14, 38), new THREE.Vector2(38, 42), new THREE.Vector2(38, 30), new THREE.Vector2(30, 2), new THREE.Vector2(30, 54), new THREE.Vector2(54, 22), new THREE.Vector2(54, 52), new THREE.Vector2(52, 20),
        new THREE.Vector2(0, 28), new THREE.Vector2(28, 52), new THREE.Vector2(28, 36), new THREE.Vector2(36, 12), new THREE.Vector2(36, 40), new THREE.Vector2(40, 16), new THREE.Vector2(16, 6), new THREE.Vector2(6, 17), new THREE.Vector2(17, 41), new THREE.Vector2(41, 37), new THREE.Vector2(56, 24),
        new THREE.Vector2(24, 10), new THREE.Vector2(10, 26), new THREE.Vector2(10, 11), new THREE.Vector2(11, 27), new THREE.Vector2(27, 51), new THREE.Vector2(40, 32), new THREE.Vector2(16, 48), new THREE.Vector2(48, 24), new THREE.Vector2(48, 49), new THREE.Vector2(49, 25), new THREE.Vector2(25, 11),
        new THREE.Vector2(25, 57), new THREE.Vector2(57, 59), new THREE.Vector2(59, 27), new THREE.Vector2(59, 35), new THREE.Vector2(35, 43), new THREE.Vector2(17, 49), new THREE.Vector2(41, 33), new THREE.Vector2(33, 1), new THREE.Vector2(1, 29), new THREE.Vector2(35, 3), new THREE.Vector2(3, 31),
        new THREE.Vector2(1, 3), new THREE.Vector2(33, 57)
    ];

    return [
        truncated_icosahedron_vertices_coordinates,
        edges_between_vertices,
        sphere_radius,
        cylinder_radius,
        sphere_rows,
        sphere_cols,
        cylinder_rows,
        cylinder_cols,
        material
    ]
}

function generate_circle(segments, radius, normal) {
    if ( radius === undefined ) radius  = 1;
    if ( normal === undefined ) normal = new THREE.Vector3(0, 0, 1);

    var angle = (2 * Math.PI) / segments;
    var vertex_buffer = new Array(segments);
    for (var i = 0; i < segments; ++i) {
        vertex_buffer[i] = new THREE.Vector3(radius*Math.cos(i*angle), radius*Math.sin(i*angle), 0);
    }

    var rotation_matrix = get_rotation_matrix(new THREE.Vector3(0, 0, 1), normal);

    for (var vertex of vertex_buffer) {
        vertex.applyMatrix4(rotation_matrix);
    }

    return vertex_buffer;
}

function generate_cylinder_geometry(first_vertex_coordinate, second_vertex_coordinate,
                               cylinder_rows, cylinder_cols, cylinder_radius, sphere_radius) {

    var distance_between_vertices = first_vertex_coordinate.clone().sub(second_vertex_coordinate).length();
    var normalized_vector = second_vertex_coordinate.clone().sub(first_vertex_coordinate).divideScalar(distance_between_vertices);
    var cylinder_length = distance_between_vertices - 1.75 * sphere_radius;
    var cylinder = new THREE.CylinderGeometry(cylinder_radius, cylinder_radius, cylinder_length,
                                                    cylinder_rows, cylinder_cols, true);
    var rotation_matrix = get_rotation_matrix(new THREE.Vector3(0, 1, 0), normalized_vector);
    cylinder.applyMatrix(rotation_matrix);
    var middle_point = first_vertex_coordinate.clone().add(second_vertex_coordinate).divideScalar(2);
    cylinder.applyMatrix(new THREE.Matrix4().makeTranslation(middle_point.x, middle_point.y, middle_point.z));
    return cylinder
}

function torus_data(rows=30, cols=15) {
    if ( rows === undefined ) rows  = 30;
    if ( cols === undefined ) cols = 15;

    var sphere_radius = 10;
    var sphere_rows = 10;
    var sphere_cols = 10;
    var cylinder_radius = 0.25 * sphere_radius;
    var cylinder_rows = 10;
    var cylinder_cols = 10;
    var material = new THREE.MeshPhongMaterial();
    material.color = new THREE.Color(0.1, 0.6, 0.1);

    var torus_vertices_coordinates_generator = function* () {
        var torus_global_radius = 8;
        var angle_local = Math.PI / cols;
        var angle_global = 2 * Math.PI / rows;
        var vector_from_base_to_radius = new THREE.Vector3(0, torus_global_radius, 0);
        var vertex = undefined;

        // Prepare normal circle
        var circle_vertices_coordinates = generate_circle(cols, 2, new THREE.Vector3(1, 0, 0));
        for (vertex of circle_vertices_coordinates) {
            vertex.add(vector_from_base_to_radius);
        }

        // Prepare shifted circle
        var shifted_circle_vertices_coordinates = generate_circle(cols, 2, new THREE.Vector3(1, 0, 0));
        var shift_rotate_matrix = new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(1, 0, 0), angle_local);
        var global_rotation_matrix = new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(0, 0, 1), angle_global/3);
        for (vertex of shifted_circle_vertices_coordinates) {
            vertex.applyMatrix4(shift_rotate_matrix);
            vertex.add(vector_from_base_to_radius);
            vertex.applyMatrix4(global_rotation_matrix);
        }

        //
        var rotation_matrix;
        for (var row = 0; row < rows; ++row) {
            if (row % 2 == 0) {
                rotation_matrix = new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(0, 0, 1), row*angle_global);
                for (vertex of circle_vertices_coordinates) {
                    yield vertex.clone().applyMatrix4(rotation_matrix);
                }
                for (vertex of shifted_circle_vertices_coordinates) {
                    yield vertex.clone().applyMatrix4(rotation_matrix);
                }
            } else {
                rotation_matrix = new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(0, 0, 1), row*angle_global - angle_global/3);
                for (vertex of shifted_circle_vertices_coordinates) {
                    yield vertex.clone().applyMatrix4(rotation_matrix);
                }

                rotation_matrix = new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(0, 0, 1), row*angle_global + angle_global/3);
                for (vertex of circle_vertices_coordinates) {
                    yield vertex.clone().applyMatrix4(rotation_matrix);
                }
            }
        }
    };

    var edges_between_vertices = function* () {
        var first_index;
        var second_index;
        var third_index;
        var fourth_index;

        for (var i = 0; i < 2*rows; i+=2) {
            for (var j = 0; j < cols; ++j) {
                first_index = i * cols + j;
                second_index = ((i + 1) * cols) + j;

                if (i % 4 == 0) {
                    third_index = ((i + 1) * cols) + Math.abs((j - 1) % cols);
                    if (j == 0) {
                        third_index +=  cols - 2;
                    }
                } else {
                    third_index = ((i + 1) * cols) + Math.abs((j + 1) % cols);
                }
                fourth_index = Math.abs(((i - 1) * cols) % (2 * rows * cols) + j);
                if (i == 0) {
                    fourth_index = 2 * rows * cols - fourth_index;
                }

                yield new THREE.Vector2(first_index, second_index);
                yield new THREE.Vector2(first_index, third_index);
                yield new THREE.Vector2(first_index, fourth_index);
            }
        }
    };

    var vertex;
    return [[for (vertex of torus_vertices_coordinates_generator()) vertex],
            edges_between_vertices(),
            sphere_radius,
            cylinder_radius,
            sphere_rows,
            sphere_cols,
            cylinder_rows,
            cylinder_cols,
            material]
}

function addMatrix(matrix1, matrix2) {
    if (matrix1.elements.length != matrix2.elements.length) {
        throw new Error( "Matrices dimensions don't match." );
    }

    var matrix_length = matrix1.elements.length;
    var result_matrix;
    if (matrix_length === 9) {
        result_matrix = new THREE.Matrix3();
    } else if (matrix_length === 16) {
        result_matrix = new THREE.Matrix4();
    } else {
        throw new Error( "Matrices with size " + matrix_length + " are not supported." );
    }

    var matrix1_elements = matrix1.elements;
    var matrix2_elements = matrix2.elements;
    var result_matrix_elements = result_matrix.elements;
    for (var i = 0; i < matrix_length; ++i) {
        result_matrix_elements.set(i, matrix1_elements.get(i) + matrix2_elements.get(i));
    }
    return result_matrix;
}

function addMatrixInplace(matrix1, matrix2) {
    if (matrix1.elements.length != matrix2.elements.length) {
        throw new Error( "Matrices dimensions don't match. "
                         + matrix1.elements.length
                         + " vs "
                         + matrix2.elements.length );
    }

    var matrix_length = matrix1.elements.length;
    if (matrix_length !== 9 && matrix_length !== 16) {
        throw new Error( "Matrices with size " + matrix_length + " are not supported." );
    }

    var matrix1_elements = matrix1.elements;
    var matrix2_elements = matrix2.elements;
    for (var i = 0; i < matrix_length; ++i) {
        matrix1_elements[i] = matrix1_elements[i] + matrix2_elements[i];
    }
    return matrix1_elements;
}

function dotMatrix(matrix1, matrix2) {
    if (matrix1.elements.length != matrix2.elements.length) {
        throw new Error( "Matrices dimensions don't match." );
    }

    var matrix_length = matrix1.elements.length;
    var result_matrix;
    if (matrix_length === 9) {
        result_matrix = new THREE.Matrix3();
    } else if (matrix_length === 16) {
        result_matrix = new THREE.Matrix4();
    } else {
        throw new Error( "Matrices with size " + matrix_length + " are not supported." );
    }

    var matrix_dimension = Math.sqrt(matrix_length);
    var matrix1_elements = matrix1.elements;
    var matrix2_elements = matrix2.elements;
    var result_matrix_elements = result_matrix.elements;
    var new_cell_value;
    for (var i = 0; i < matrix_dimension; ++i) {
        for (var j = 0; j < matrix_dimension; ++j) {
            new_cell_value = 0;
            for (var k = 0; k < matrix_dimension; ++k) {
                new_cell_value += matrix1_elements[i*matrix_dimension + k] * matrix2_elements[k*matrix_dimension + j];
            }
            result_matrix_elements[i*matrix_dimension + j] = new_cell_value;
        }
    }

    return result_matrix;
}

function matrix4toMatrix3(matrix4) {
    var matrix3 = new THREE.Matrix3();
	var matrix4_elements = matrix4.elements;

    matrix3.fromArray([
        matrix4_elements[0], matrix4_elements[1], matrix4_elements[2],
        matrix4_elements[4], matrix4_elements[5], matrix4_elements[6],
        matrix4_elements[8], matrix4_elements[9], matrix4_elements[10]
    ]);

    return matrix3
}

function matrix3toMatrix4(matrix3) {
    var matrix4 = new THREE.Matrix4();
	var matrix3_elements = matrix3.elements;

    matrix4.fromArray([
        matrix3_elements[0], matrix3_elements[1], matrix3_elements[2], 0,
        matrix3_elements[3], matrix3_elements[4], matrix3_elements[5], 0,
        matrix3_elements[6], matrix3_elements[7], matrix3_elements[8], 0,
						  0, 				   0,					0, 1
    ]);

    return matrix4
}

function get_rotation_matrix(start_point, end_point) {
    if (start_point.equals(end_point)) {
        return new THREE.Matrix4();
    }

    var x = start_point.clone().cross(end_point).normalize();
    var start_point_length = start_point.length();
    var end_point_length = end_point.length();
    var cos_theta = start_point.dot(end_point) / (start_point_length * end_point_length);
    var theta = Math.acos(cos_theta);
    var A = new THREE.Matrix3().fromArray(
        [   0, -x.z,  x.y,
          x.z,    0, -x.x,
         -x.y,  x.x,    0 ]);
    var rotate_matrix = new THREE.Matrix3();
    addMatrixInplace(rotate_matrix, A.clone().multiplyScalar(Math.sin(theta)));
    addMatrixInplace(rotate_matrix, dotMatrix(A, A).multiplyScalar(1 - cos_theta));

    if (x.z != 0) {
        rotate_matrix = dotMatrix(rotate_matrix, matrix4toMatrix3(new THREE.Matrix4().makeRotationZ(Math.PI)));
    }
    if (x.y != 0) {
        rotate_matrix = dotMatrix(rotate_matrix, matrix4toMatrix3(new THREE.Matrix4().makeRotationY(Math.PI)));
    }
    rotate_matrix = dotMatrix(rotate_matrix, matrix4toMatrix3(new THREE.Matrix4().makeRotationX(Math.PI)));

    return matrix3toMatrix4(rotate_matrix);
}

function generate_skeleton(vertices_coordinates_generator, edges_between_vertices_generator, sphere_radius,
                           cylinder_radius, sphere_rows, sphere_cols, cylinder_rows, cylinder_cols, material)
{
    var final_object_geometry = new THREE.Geometry();
    var object_geometry;
    var current_vertex;
    for (var vertex of vertices_coordinates_generator) {
        current_vertex = vertex.clone().multiplyScalar(100);
        object_geometry = new THREE.SphereGeometry(sphere_radius, sphere_rows, sphere_cols);
        object_geometry.applyMatrix(new THREE.Matrix4().makeTranslation(current_vertex.x,
                                                                        current_vertex.y,
                                                                        current_vertex.z));
        final_object_geometry.merge(object_geometry);
    }

    var first_vertex_coordinate;
    var second_vertex_coordinate;
    for (var edge of edges_between_vertices_generator) {
        first_vertex_coordinate = vertices_coordinates_generator[edge.x].clone().multiplyScalar(100);
        second_vertex_coordinate = vertices_coordinates_generator[edge.y].clone().multiplyScalar(100);
        object_geometry = generate_cylinder_geometry(first_vertex_coordinate, second_vertex_coordinate,
                                                     cylinder_rows, cylinder_cols, cylinder_radius, sphere_radius);
        final_object_geometry.merge(object_geometry);
    }

    return new THREE.Mesh( final_object_geometry, material );
}

function init() {
    scene = new THREE.Scene();
    //
    var mesh = generate_skeleton.apply(this, truncated_icosahedron_vertices_coordinates_data());
    scene.add( mesh );
    //
    var mesh = generate_skeleton.apply(this, torus_data(30, 15));
    scene.add( mesh );
	//
	var axisHelper = new THREE.AxisHelper( 100 );
	scene.add( axisHelper );
    //
    var light = new THREE.DirectionalLight( 0xffffff, 1);
    light.position.set( 50, 50, 50 ).normalize();
    scene.add( light );
    //
    renderer = new THREE.WebGLRenderer();
    renderer.setPixelRatio( window.devicePixelRatio );
    renderer.setSize( window.innerWidth, window.innerHeight );
    document.body.appendChild( renderer.domElement );
	//
	camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 1, 3000 );
    camera.position.x = 400;
    camera.position.y = 400;
    camera.position.z = 1000;
    camera.up = new THREE.Vector3(0,1,0);
    camera.lookAt(new THREE.Vector3(0,0,0));
	cameraControls = new THREE.OrbitControls(camera, renderer.domElement);
	cameraControls.target.set( 0, 0, 0);
	cameraControls.maxDistance = 3000;
	cameraControls.minDistance = 10;
	cameraControls.update();
    //
    window.addEventListener( 'resize', onWindowResize, false );
}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
}

function animate() {
    requestAnimationFrame( animate );
    scene.children[0].rotation.x += 0.005;
    scene.children[0].rotation.y += 0.01;
    scene.children[1].rotation.x -= 0.005;
    scene.children[1].rotation.y -= 0.01;
	cameraControls.update();
    renderer.render( scene, camera );
}