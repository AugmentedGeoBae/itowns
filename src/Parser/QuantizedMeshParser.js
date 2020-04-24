/*
*Quantized Mesh is Triangulated Irregular Network (TIN) mesh format
*
*https://github.com/CesiumGS/quantized-mesh
*
*/



import * as THREE from './build/three.module.js';



function highwaterDecode(indices) {

    //var arr = [];

	//
	// use same array as input
	//

	console.log("indices in ", indices);

    let highest = 0;
    for (let i = 0; i < indices.length; ++i) {
        let code = indices[i];
		//console.log("code ", code);
        indices[i] = highest - code;
		//console.log("arr[i] ", arr[i]);
		//arr.push(highest - code);
        if (code === 0) {
            ++highest;
        }
    }

	console.log("indices out ", indices);
    return indices;
}

function getEdgeIndices(data, indPos, ind_vCount, byteCount) {

	let edgeInds;

	if (byteCount === 4) {
		edgeInds = getUint32Array(data.buffer, indPos, ind_vCount);
		//indArr = getUint32Array(data.buffer, indPos, ind_vCount)
		//edgeInds = highwaterDecode(indArr);
	} else {
		edgeInds = getUint16Array(data.buffer, indPos, ind_vCount);
		//indArr = getUint16Array(data.buffer, indPos, ind_vCount)
		//edgeInds = highwaterDecode(indArr, returnArr);
	}

	return edgeInds;

}

function getEdges(data, edgePos, n_vertices) {

  let edges = {};

  let n_bytes;

  if (n_vertices > 65536) {
	n_bytes = 4;
  } else {
	n_bytes = 2;
  }

  const west_vCount = data.getUint32(edgePos, true);
  edgePos += Uint32Array.BYTES_PER_ELEMENT;
  const w_indices = getEdgeIndices(data, edgePos, west_vCount, n_bytes);
  edgePos += west_vCount * n_bytes;

  const south_vCount = data.getUint32(edgePos, true);
  edgePos += Uint32Array.BYTES_PER_ELEMENT;
  const s_indices = getEdgeIndices(data, edgePos, south_vCount, n_bytes);
  edgePos += south_vCount * n_bytes;

  const east_vCount = data.getUint32(edgePos, true);
  edgePos += Uint32Array.BYTES_PER_ELEMENT;
  const e_indices = getEdgeIndices(data, edgePos, east_vCount, n_bytes);
  edgePos += east_vCount * n_bytes;

  const north_vCount = data.getUint32(edgePos, true);
  edgePos += Uint32Array.BYTES_PER_ELEMENT;
  const n_indices = getEdgeIndices(data, edgePos, north_vCount, n_bytes);
  edgePos += north_vCount * n_bytes;

  //let indicesDecoded = nOctDecode(edgeIndices);

  console.log("n_indices ", n_indices, n_indices[2]);

  edges.west = w_indices;
  edges.south = s_indices;
  edges.east = e_indices;
  edges.north = n_indices;

  console.log("edges ", edges, edgePos);

  return {edges, edgePos};

}

function signNotZero (vec) {
  return new THREE.Vector2( vec.x >= 0 ? 1 : -1, vec.y >= 0 ? 1 : -1 )
}

function nOctDecode (encodedVec) {
  let decodedVec = encodedVec.divideScalar(255).multiplyScalar(2).subScalar(1);

  decodedVec = new THREE.Vector3( decodedVec.x, decodedVec.y, 1 - Math.abs(decodedVec.x) - Math.abs(decodedVec.y) )

  if (decodedVec.z < 0) {
    const xy = new THREE.Vector2(decodedVector.x, decodedVector.y)
    const xyAbs = xy.distanceTo(new THREE.Vector2(0, 0))
    const xySign = signNotZero(xy)
    const decodedXy = xySign.multiplyScalar(1 - xyAbs)

    decodedVec.set(decodedXy.x, decodedXy.y, decodedVec.z)
  }

  return decodedVec.normalize()
}

const getNormals = (buff, data, startPos, normalsLen) => {

  //normals are x,y and 2 byte encoded
  let normalArr = new Uint8Array(buff, starPos, normalsLen);

  const view = new DataView(vertexNormalsBuffer)
  const elementsPerEncodedNormal = 2
  const elementsPerNormal = 3
  const vertexNormalsAttributeArray = new Float32Array(vertexData.length)

  for (let position = 0, i = 0; position < vertexNormalsBuffer.byteLength; position += Uint8Array.BYTES_PER_ELEMENT * elementsPerEncodedNormal, i++) {
    const decodedNormal = decodeOct(new THREE.Vector2(
      view.getUint8(position, true),
      view.getUint8(position + Uint8Array.BYTES_PER_ELEMENT, true)
    ))

    vertexNormalsAttributeArray[i * elementsPerNormal] = decodedNormal.x
    vertexNormalsAttributeArray[i * elementsPerNormal + 1] = decodedNormal.y
    vertexNormalsAttributeArray[i * elementsPerNormal + 2] = decodedNormal.z
  }

  //return new THREE.BufferAttribute(vertexNormalsAttributeArray, elementsPerNormal, true)

  let vNormals = vertexNormalsAttributeArray;

  return vNormals;
}

const getWaterMask = (data, startPos, maskLen) => {

  let waterMask = data.buffer.slice(startPos, startPos + maskLen);

  let type;

  if (waterMask.length() > 1) {
    type = "mixed";
  } else {
    if (data.getUint8(waterMask[0]) === 0) {
      type = "land";
    } else {
      type = "water";
    }
  }

  return waterMask, type;
}

const getMetadata = (data, startPos, metadataLen) => {

  let metadataPos = startPos;
  let jsonLen = data.getUint8(startpos, true);
  metadataPos += Uint8Array.BYTES_PER_ELEMENT;

  let metadataJson = JSON.parse(data.buffer.slice(metadataPos, metadataLen));

  return metadataJson;

}

const zigZagDecode = (value) => {
    return (value >> 1) ^ (-(value & 1));
}

const getUint32Array = (data, startPos, count) => {
    return new Uint32Array(data.slice(startPos, startPos + 4*count));
}

const getUint16Array = (data, startPos, count) => {
	//console.log("getUint16Array ", data, startPos, count);
    return new Uint16Array(data.slice(startPos, startPos + Uint16Array.BYTES_PER_ELEMENT*count));
}

const getFloat32Array = (data, startPos, count) => {
    return new Float32Array(data.slice(startPos, startPos + 4*count));
}

const getHeader = (data, byteOffset) => {

  return {
      bytes: data.byteLength,
      centerX : data.getFloat64(byteOffset, true), //getFloat64(data, byteOffset),
      centerY : data.getFloat64(byteOffset + 8, true),
      centerZ : data.getFloat64(byteOffset + 16, true),
      minimumHeight : data.getFloat32(byteOffset + 24, true),
      maximumHeight : data.getFloat32(byteOffset + 28, true),
      boundingSphereCenterX : data.getFloat64(byteOffset + 32, true),
      boundingSphereCenterY : data.getFloat64(byteOffset + 40, true),
      boundingSphereCenterZ : data.getFloat64(byteOffset + 48, true),
      boundingSphereRadius : data.getFloat64(byteOffset + 56, true),
      horizonOcclusionPointX : data.getFloat64(byteOffset + 64, true),
      horizonOcclusionPointY : data.getFloat64(byteOffset + 72, true),
      horizonOcclusionPointZ : data.getFloat64(byteOffset + 80, true)
  }

}

function getVertices(data, startPos, vCount) {

	let vPos = startPos;

	const uArray = getUint16Array(data.buffer, vPos, vCount);
    vPos += vCount * Uint16Array.BYTES_PER_ELEMENT;

	console.log("uArray ", uArray.byteLength);

    const vArray = getUint16Array(data.buffer, vPos, vCount);
    vPos += vCount * Uint16Array.BYTES_PER_ELEMENT;

    const heightArray = getUint16Array(data.buffer, vPos, vCount);
    vPos += vCount * Uint16Array.BYTES_PER_ELEMENT;

	console.log("byteOffset vertexCount*2*3 ", vPos, data.byteLength); //correct

    let i;
    let u = 0;
    let v = 0;
    let height = 0;

    for (i = 0; i < uArray.length; ++i) {
        u += zigZagDecode(uArray[i]);
        v += zigZagDecode(vArray[i]);
        height += zigZagDecode(heightArray[i]);

        uArray[i] = u;
        vArray[i] = v;
        heightArray[i] = height;
    }

	// padding for byte alignment?
    if (vPos % 2 !== 0) {
        vPos += (2 - (vPos % 2));
    }

    let vertices = new Uint16Array(vCount*3);

    for (let j = 0; j < uArray.length; j++) {
        //vertices.push(new THREE.Vector3(uArray[i]/100, vArray[i]/100, heightArray[i]/200));
        vertices[j] = (uArray[j]);
        vertices[j + vCount] = (vArray[j]);
        vertices[j + vCount*2] = (heightArray[j]);
    }

	console.log("vertices ", vertices);

    return {vertices, vPos};
};

function getFaces(uArray, vArray, heightArray, indexArray) {
    var faces = [];
    for (var i = 0; i < indexArray.length; i+=3) {
        faces.push(new THREE.Face3(indexArray[i], indexArray[i+1], indexArray[i+2]));
    }
    return faces;
};

function getGeometry(indexArray, vertices) {
  let geometry = new THREE.BufferGeometry();

  console.log("vertices at getGeometry ", vertices);

  let posArr = new Float32Array(vertices.length);

  for (let k = 0; k < vertices.length; k++) {

	posArr[k * 3] = vertices[k];
	posArr[k * 3 + 1] = vertices[k + (vertices.length/3)];
	posArr[k * 3 + 2] = vertices[k + (2* vertices.length/3)];

  }

  geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( posArr, 3 ) );

  geometry.setIndex( new THREE.BufferAttribute(indexArray, 1) );
  console.log("geometry ", geometry);

  return geometry;
};

function getExtensions(data, extPos) {

  let exts = [];
  //let byteOffset1 = byteOffset;

  while (extPos < data.byteLength) {

    var extId = data.getUint8(extPos, true); //1 byte (8 bits)
    extPos += Uint8Array.BYTES_PER_ELEMENT;
    var extLength = data.getUint32(extPos, true);
    extPos += Uint32Array.BYTES_PER_ELEMENT;

    extPos += extLength;

    if (extId == 1) {

      exts.push( getNormals(data, extPos, extLength) );

    } else if (extId == 2) {

      exts.push( getWaterMask(data, extPos, extLength) );

    } else if (extId == 4) {

	  exts.push( getMetadata(data, extPos, extLength) );

	} else {

	  console.warn(`ID ${extId} not recognised, only 1[= vertex normals] and 2[= water mask] are accepted`);

    }

  }

  return { exts, extPos };

}

export default {
  /** Parse quantized mesh buffer and extract THREE.Scene and batch table
   * @param {ArrayBuffer} buffer - the quantized mesh buffer.
   * @param {Object} options - additional properties.
   * @param {string=} [options.quantizedMeshUpAxis='Y'] - mesh up axis.
   * @param {string} options.urlBase - the base url of the b3dm file (used to fetch textures for the embedded glTF model).
   * @param {boolean=} [options.doNotPatchMaterial='false'] - disable patching material with logarithmic depth buffer support.
   * @param {float} [options.opacity=1.0] - the b3dm opacity.
   * @param {boolean|Material=} [options.overrideMaterials='false'] - override b3dm's embedded glTF materials. If overrideMaterials is a three.js material, it will be the material used to override.
   * @return {Promise} - a promise that resolves with an object containig a THREE.Scene (gltf) and a batch table (batchTable).
   *
   */

  //get path url then convert to arrayBuffer : pathToFile.arrayBuffer()

  parse: function parse(buffer, options) {

    const view = new DataView(buffer);

    let quantizedMeshComponents = {};

    let byteOffset = 0;

    quantizedMeshComponents.header = getHeader(view, byteOffset);

    byteOffset += 88;

	console.log("byteOffset 88 ", byteOffset);

    const vertexCount = view.getUint32(88, true);
    byteOffset += Uint32Array.BYTES_PER_ELEMENT;

	console.log("byteOffset 92 ", byteOffset);
	console.log("vertexCount ", vertexCount);
	console.log("view length ", view.buffer.byteLength, view.byteLength);

	const {vertices, vPos: verticesEndPos} = getVertices(view, byteOffset, vertexCount);

	byteOffset = verticesEndPos;

	console.log("verticesEndPos, vertices ", byteOffset, vertices);

    const triangleCount = view.getUint32(byteOffset, true);
    byteOffset += Uint32Array.BYTES_PER_ELEMENT;

	const indicesCount = triangleCount *3;
    //To enforce proper byte alignment,
    //padding is added before the IndexData to ensure 2 byte alignment for IndexData16 and 4 byte alignment for IndexData32.

    let indices;

    console.log("vertex count ", vertexCount);
    console.log("triangleCount ", triangleCount);

    if (vertexCount > 65536) {

      indices = getUint32Array(buffer, byteOffset, indicesCount);

      byteOffset += indicesCount * Uint32Array.BYTES_PER_ELEMENT;

    } else {

      indices = getUint16Array(buffer, byteOffset, indicesCount);

      byteOffset += indicesCount * Uint16Array.BYTES_PER_ELEMENT;

    }

	console.log("after indices pos ", byteOffset, indices.byteLength);

    const indexArray = highwaterDecode(indices);

/////////////////////////////////

    //why does it calculate faces?
    //const faces = getFaces(uArray, vArray, heightArray, indexArray);

    quantizedMeshComponents.geometry = getGeometry(indexArray, vertices);

/////////////////////////////////

    // triangleCount * 3 = vertex count?
	//let quantizedMeshEdges, edgeEndPos;
    const quantizedMeshEdges = getEdges(view, byteOffset, vertexCount);

	console.log("quantizedMeshEdges ", quantizedMeshEdges);

	quantizedMeshComponents.edges = quantizedMeshEdges.edges;//quantizedMeshEdges;

	//let edgeEndPos = bcde;

	console.log("edgeEndPos ", quantizedMeshEdges.edgePos);

	byteOffset = quantizedMeshEdges.edgePos; //edgeEndPos;

	//let extEndPos;
	//{ quantizedMeshComponents.extensions, extEndPos }
    const quantizedMeshExts = getExtensions(view, byteOffset);

	quantizedMeshComponents.extensions = quantizedMeshExts.exts;

    byteOffset = quantizedMeshExts.extPos;

	console.log("byteOffset end ", byteOffset, view.byteLength);

    if (byteOffset != view.byteLength) {

      console.warn("byteOffset != view.byteLength");

    }

    //return Promise.resolve(quantizedMeshTile_bufferGeometry);

    console.log(quantizedMeshComponents);

    return Promise.resolve(quantizedMeshComponents);

	//return quantizedMeshComponents;

    //It should return a buffer geometry as part of an object containing vertexes, indices, edges, extensions

  }
};
