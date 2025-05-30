/*
 *  Copyright (c), 2017, Adrien Devresse <adrien.devresse@epfl.ch>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 *
 */
#pragma once

#include "../H5Easy.hpp"

namespace H5Easy {

namespace detail {

// Generate error-stream and return "Exception" (not yet thrown).
inline Exception error(const File& file, const std::string& path, const std::string& message) {
    std::ostringstream ss;
    ss << message << std::endl
       << "Path: " << path << std::endl
       << "Filename: " << file.getName() << std::endl;
    return Exception(ss.str());
}

// Generate specific dump error
inline Exception dump_error(File& file, const std::string& path) {
    if (file.getObjectType(path) == ObjectType::Dataset) {
        return error(file,
                     path,
                     "H5Easy: Dataset already exists, dump with H5Easy::DumpMode::Overwrite "
                     "to overwrite (with an array of the same shape).");
    } else {
        return error(
            file,
            path,
            "H5Easy: path exists, but does not correspond to a Dataset. Dump not possible.");
    }
}

// get a opened DataSet: nd-array
template <class T>
inline DataSet initDataset(File& file,
                           const std::string& path,
                           const std::vector<size_t>& shape,
                           const DumpOptions& options) {
    if (!file.exist(path)) {
        if (!options.compress() && !options.isChunked()) {
            return file.createDataSet<T>(path, DataSpace(shape), {}, {}, true);
        } else {
            std::vector<hsize_t> chunks(shape.begin(), shape.end());
            if (options.isChunked()) {
                chunks = options.getChunkSize();
                if (chunks.size() != shape.size()) {
                    throw error(file, path, "H5Easy::dump: Incorrect rank ChunkSize");
                }
            }
            DataSetCreateProps props;
            props.add(Chunking(chunks));
            if (options.compress()) {
                props.add(Shuffle());
                props.add(Deflate(options.getCompressionLevel()));
            }
            return file.createDataSet<T>(path, DataSpace(shape), props, {}, true);
        }
    } else if (options.overwrite() && file.getObjectType(path) == ObjectType::Dataset) {
        DataSet dataset = file.getDataSet(path);
        if (dataset.getDimensions() != shape) {
            throw error(file, path, "H5Easy::dump: Inconsistent dimensions");
        }
        return dataset;
    }
    throw dump_error(file, path);
}

// get a opened DataSet: scalar
template <class T>
inline DataSet initScalarDataset(File& file,
                                 const std::string& path,
                                 const T& data,
                                 const DumpOptions& options) {
    if (!file.exist(path)) {
        return file.createDataSet<T>(path, DataSpace::From(data), {}, {}, true);
    } else if (options.overwrite() && file.getObjectType(path) == ObjectType::Dataset) {
        DataSet dataset = file.getDataSet(path);
        if (dataset.getElementCount() != 1) {
            throw error(file, path, "H5Easy::dump: Existing field not a scalar");
        }
        return dataset;
    }
    throw dump_error(file, path);
}

template <class File, class F>
auto apply_attr_func_impl(File& file, const std::string& path, F f) {
    auto type = file.getObjectType(path);
    if (type == ObjectType::Group) {
        auto group = file.getGroup(path);
        return f(group);
    } else if (type == ObjectType::Dataset) {
        auto dataset = file.getDataSet(path);
        return f(dataset);
    } else {
        throw error(file, path, "path is not the root, a group or a dataset.");
    }
}

template <class F>
auto apply_attr_func(const H5Easy::File& file, const std::string& path, F f) {
    return apply_attr_func_impl(file, path, f);
}

template <class F>
auto apply_attr_func(H5Easy::File& file, const std::string& path, F f) {
    return apply_attr_func_impl(file, path, f);
}

// get a opened Attribute: nd-array
template <class T>
inline Attribute initAttribute(File& file,
                               const std::string& path,
                               const std::string& key,
                               const std::vector<size_t>& shape,
                               const DumpOptions& options) {
    auto get_attribute = [&](auto& obj) {
        if (!obj.hasAttribute(key)) {
            return obj.template createAttribute<T>(key, DataSpace(shape));
        } else if (options.overwrite()) {
            Attribute attribute = obj.getAttribute(key);
            DataSpace dataspace = attribute.getSpace();
            if (dataspace.getDimensions() != shape) {
                throw error(file, path, "H5Easy::dumpAttribute: Inconsistent dimensions");
            }
            return attribute;
        }
        throw error(file,
                    path,
                    "H5Easy: Attribute exists, overwrite with H5Easy::DumpMode::Overwrite.");
    };

    if (!file.exist(path)) {
        throw error(file, path, "H5Easy::dumpAttribute: path does not exist");
    }

    return apply_attr_func(file, path, get_attribute);
}

// get a opened Attribute: scalar
template <class T>
inline Attribute initScalarAttribute(File& file,
                                     const std::string& path,
                                     const std::string& key,
                                     const T& data,
                                     const DumpOptions& options) {
    auto get_attribute = [&](auto& obj) {
        if (!obj.hasAttribute(key)) {
            return obj.template createAttribute<T>(key, DataSpace::From(data));
        } else if (options.overwrite()) {
            Attribute attribute = obj.getAttribute(key);
            DataSpace dataspace = attribute.getSpace();
            if (dataspace.getElementCount() != 1) {
                throw error(file, path, "H5Easy::dumpAttribute: Existing field not a scalar");
            }
            return attribute;
        }
        throw error(file,
                    path,
                    "H5Easy: Attribute exists, overwrite with H5Easy::DumpMode::Overwrite.");
    };

    if (!file.exist(path)) {
        throw error(file, path, "H5Easy::dumpAttribute: path does not exist");
    }

    apply_attr_func(file, path, get_attribute);
}

}  // namespace detail
}  // namespace H5Easy
